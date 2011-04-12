/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/luminaire.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <set>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

MTS_NAMESPACE_BEGIN

/**
 * Height Span Map loader
 */
class HeightSpanMap : public Shape {
public:

    /**
     * Represents an height element of within the height span of a cell.
     */
    struct HeightSample
    {
        /* The type of this heigt span element:
         * Normal = 0, Bridge = 1, Bridge border = 2
         */
        int type;
        // Tha base height of this element.
        float h;
        // The snow height on this element, relative to h.
        float dh;
        // North and east offset of snow node. Used for bridges and over hangs
        float dx, dy;
        //
        float u, v;
        // Flags to state in which direction a neighbor exists.
        unsigned char flags, user_flags;
        // indices of adjacent spans
        short neighbor_slab_indices[4];
        /* Distances to adjacent spans. Indexed after incoming directions.
         * This means: If you are interested in distance to the neighbor
         * east of this element, you need to ask for the west distance.
         */
        short neighbor_distances[4];

        int vertex_index, secondary_vertex_index;

        /**
         * Creates a new height span element. Defaults are at height 0, no
         * span height and no neighbors.
         */
        HeightSample(float _h = 0, float _dh = 0, unsigned char _f = 0)
            : type(0),h(_h), dh(_dh), dx(0), dy(0), u(0), v(0), flags(_f), user_flags(0), vertex_index(-1), secondary_vertex_index(-1)
        {
            // init neighbor slab indices to -1
            neighbor_slab_indices[0] = neighbor_slab_indices[1] = neighbor_slab_indices[2] = neighbor_slab_indices[3] = -1;
            // init neighbor distances to 1
            neighbor_distances[0] = neighbor_distances[1] = neighbor_distances[2] = neighbor_distances[3] = 1;
        }

        /**
         * Creates a new heaght span element. Defaults to neighborless.
         */
        HeightSample(float _h, int _type, float _dh, float _dx, float _dy, unsigned char _f = 0)
            : type(_type), h(_h), dh(_dh), dx(_dx), dy(_dy), u(0), v(0), flags(_f), user_flags(0), vertex_index(-1), secondary_vertex_index(-1)
        {
            // init neighbor slab indices to -1
            neighbor_slab_indices[0] = neighbor_slab_indices[1] = neighbor_slab_indices[2] = neighbor_slab_indices[3] = -1;
            // init neighbor distances to 1
            neighbor_distances[0] = neighbor_distances[1] = neighbor_distances[2] = neighbor_distances[3] = 1;
        }

        /**
         * Return absolut height of this element. Tat is base height
         * plus snow height.
         */
        float get_total_height() const {
            return h + dh;
        }

        /**
         * Gets an extendet type information by including
         * info about whether this elements snow node is offset.
         * If there is an offset, return actual type + 3,
         * otherwise return real type.
         */
        int get_extended_type() const {
            return dx*dx+dy*dy > 0 ? type + 3 : type;
        }

        /**
         * Test if this element is a border element.
         */
        bool is_border() const {
            // return true if there is at least one direction without neighbor.
            return flags != 15;
        }

        /**
         * Check if this element has a neighbor in direction i.
         */
        bool get_flag(int i) const {
            static unsigned char flag_values[] = { 1,2,4,8 };
            return (flags & flag_values[i]) != 0;
        }

        /**
         * State that there is a neihgbor in direction i.
         */
        void set_flag(int i) {
            static unsigned char flag_values[] = { 1,2,4,8 };
            flags |= flag_values[i];
        }

        bool get_user_flag(int i) const {
            static unsigned char flag_values[] = { 1,2,4,8 };
            return (user_flags & flag_values[i]) != 0;
        }
        void set_user_flag(int i) {
            static unsigned char flag_values[] = { 1,2,4,8 };
            user_flags |= flag_values[i];
        }

        /**
         * Get the index of an adjacent element within neighboring
         * cell in direction i. It is therefore the snow slab in this
         * adjacent cell.
         */
        short get_nbr_slab_idx(int i) const {
            return neighbor_slab_indices[i];
        }

        /**
         * Get dintance to an adjacent element in direction i.
         */
        short get_distance(int i) const {
            return neighbor_distances[i];
        }

        std::string toString() {
            std::stringstream ss;
            ss << "Height sample:" << std::endl
               << "\ttype: " << type << std::endl
               << "\th: " << h << std::endl
               << "\tdh: " << dh << std::endl
               << "\tdx: " << dx << " dy: " << dy <<  std::endl
               << "\tu: " << u << " v: " << v << std::endl
               << "\tflag_0: " << get_flag(0) << " flag_1: " << get_flag(1)
               << " flag_2: " << get_flag(2) << " flag_3: " << get_flag(3) << std::endl
               << "\tnsl_0: " << neighbor_slab_indices[0] << " nsl_1: " << neighbor_slab_indices[1]
               << " nsl_2: " << neighbor_slab_indices[2] << " nsl_2: " << neighbor_slab_indices[3] << std::endl
               << "\tnd_0: " << neighbor_distances[1] << " nd_1: " << neighbor_distances[1]
               << " nd_2: " << neighbor_distances[2] << " nd_3: " << neighbor_distances[3] << std::endl
               << "\tvertex-idx: " << vertex_index << " sec. vertex-idx: " << secondary_vertex_index << std::endl;

            return ss.str();
        }
    };

    /**
     * Data structure for keeping position information on the map.
     */
    struct HeightSpanIterator
    {
        float x,y;
        float dx,dy;
        int i,j,k;
        HeightSpanIterator(float _dx, float _dy, int _i, int _j, int _k = 0)
            : x((_i/*+0.5f*/)*_dx), y((_j/*+0.5f*/)*_dy), dx(_dx), dy(_dy), i(_i), j(_j), k(_k) {}
    };

public:
    /* type definitons */
    typedef boost::char_separator<char> separator;
    typedef boost::tokenizer< separator > tokenizer;
    // Represent height spans as vector collection of height samples.
	typedef std::vector<HeightSample> HeightSpanType;


	bool fetch_line(std::istream &is, std::string &line) {
		/// Fetch a line from the stream, while handling line breaks with backslashes
		if (!std::getline(is, line))
			return false;
		if (line == "")
			return true;
		int lastCharacter = line.size()-1;
		while (lastCharacter >= 0 &&
				(line[lastCharacter] == '\r' ||
					line[lastCharacter] == '\n' ||
					line[lastCharacter] == '\t' ||
					line[lastCharacter] == ' '))
			lastCharacter--;

		if (line[lastCharacter] == '\\') {
			std::string nextLine;
			fetch_line(is, nextLine);
			line = line.substr(0, lastCharacter) + nextLine;
		} else {
			line.resize(lastCharacter+1);
		}
		return true;
	}

	HeightSpanMap(const Properties &props) : Shape(props) {
		FileResolver *fResolver = Thread::getThread()->getFileResolver();
		fs::path path = fResolver->resolve(props.getString("filename"));
		m_name = path.stem();

		/* By default, any existing normals will be used for
		   rendering. If no normals are found, Mitsuba will
		   automatically generate smooth vertex normals. 
		   Setting the 'faceNormals' parameter instead forces
		   the use of face normals, which will result in a faceted
		   appearance.
		*/
		m_faceNormals = props.getBoolean("faceNormals", true);

		/* Re-center & scale all contents to move them into the
		   AABB [-1, -1, -1]x[1, 1, 1]? */
		m_recenter = props.getBoolean("recenter", false);

		/* Causes all normals to be flipped */
		m_flipNormals = props.getBoolean("flipNormals", false);

		/* Object-space -> World-space transformation */
		Transform objectToWorld = props.getTransform("toWorld", Transform());

		/* Causes all normals to be flipped */
		std::string name = props.getString("name", m_name);

		BSDF *currentMaterial = NULL;
	    bool hasNormals = false;
        bool hasTexcoords = false;
        
        /* Read in the data and build up metsuba usable geometry */
        read(path);
        generateGeometry(name, hasNormals, hasTexcoords,
            currentMaterial, objectToWorld);
    }

    /**
     * Reads in a file. It (re)initialises the object and does
     * version handling.
     */
    void read(fs::path path) {
		/* Load the geometry */
		Log(EInfo, "Loading height span map \"%s\" ..", path.leaf().c_str());
		fs::ifstream is(path);
		if (is.bad() || is.fail())
			Log(EError, "Height span map file '%s' not found!", path.file_string().c_str());

		std::string buf;

        // do a simple extension test to determine version
        std::string ext = path.extension();
        boost::to_lower(ext);
        int version;

        if (ext.compare(".hspans1") == 0)
            version = 1;
        else if (ext.compare(".hspans2") == 0)
            version = 2;
        else
            version = 3;

		std::string name = m_name, line;
        bool dimensionsFound = false;

        // clear the map
        m_heightSpans.clear();

        // eat all leading comments
		while (is.good() && !is.eof() && fetch_line(is, line)) {
            // skip comments
            if (line[0] == '#')
               continue;

			std::istringstream iss(line);

            int w, h;
            iss >> w >> h;
            resize(w, h);
            dimensionsFound = true;
            break;
        }

        // complain if dimensiorn were not yet found
        if (!dimensionsFound) {
            Log(EError, "Encountered an error while loading height span map: Expected dimensions");
        }

        /* in version 1 and 2 we need to skip 3 lines, 7 in version 3
         *  until the bounding box is found
         */
        int skips = (version == 3) ? 7 : 3;
		while (skips != 0 && is.good() && !is.eof() && fetch_line(is, line))
            --skips;

        // complain if skipping went wrong
        if (skips > 0) {
            Log(EError, "Encountered an error while loading height span map: Expected more data");
        }

        /* next a bounding box can be specified, but this
         * is not used or implemented yet. The bounding box is
         * always recalculated:
         */
        int bbCount = 4;
        std::vector<std::string> bbEntries(4);
        do {
            --bbCount;
            bbEntries[3 - bbCount] = line;
        } while (bbCount != 0 && is.good() && !is.eof() && fetch_line(is, line));

        // complain if we got not all bounding box entries
        if (bbCount > 0) {
            Log(EError, "Expected more bounding box data");
        }
        // read in bounding box data
        Vector o = parseLineWithVector(bbEntries[0]);
        Vector x = parseLineWithVector(bbEntries[1]) - o;
        Vector y = parseLineWithVector(bbEntries[2]) - o;
        Vector z = parseLineWithVector(bbEntries[3]) - o;
        /* use bounding box as matrix transformation:
         * translation and scaling
         */
        Matrix4x4 trafo(
            x.x, y.x, z.x, o.x,
            x.y, y.y, z.y, o.y,
            x.z, y.z, z.z, o.z,
            0, 0, 0, 1
        );

        m_fileToObject = Transform(trafo);

        // do version specific parsing
        if (version == 1)
            readDataV1(is);
        else if (version == 2)
            readDataV2(is);
        else
            readDataV3(is);
	}

    Vector parseLineWithVector(const std::string& line) const
    {
        // create zero vector (0,0,0,1) as default
        Vector v(0.0, 0.0, 0.0);
        // tokenize the string, seperate it by ":" into "tokens"
        separator sep1(":");
        // tokenize the line
        tokenizer tokens(line, sep1);
        // access the tokens through a vector
        std::vector<std::string> vec1;
        vec1.assign(tokens.begin(),tokens.end());
        // check if we have more or less than two tokens
        if (vec1.size() != 2) {
            // complain about not having found seperator and return default vector
            Log(EWarn, "Error in parsing line with vector: no \":\" found");
            return v;
        }
        // set weightspace and tokenize into "nums"
        separator sep2(" \t(,)");
        tokenizer numsTokens(vec1[1], sep2);
        // access the tokens through a vector
        std::vector<std::string> nums;
        nums.assign(numsTokens.begin(),numsTokens.end());
        // make sure we have three new tokens
        if (nums.size() != 3) {
            // no, we don't, complain about it and return default vector
            Log(EWarn, "Error in parsing line with vector: no 3 numbers found");
            return v;
        }
        // try to read every token as number component, if not possible use default
        for (int i = 0; i < 3; ++i) {
            try {
                v[i] = boost::lexical_cast<double>(nums[i]);
            } catch(boost::bad_lexical_cast &) {
                std::stringstream str;
                str << "Error in parsing line with vector: number " << i << " invalid";
                Log(EWarn, str.str().c_str());
                continue;
            }
        }
        // return vector
        return v;
    }

    void readDataV1(fs::ifstream& is) {
            Log(EError,  "Hspan 1 format not yet implemented!");
    }

    void readDataV2(fs::ifstream& is) {
        separator sep("\t ");
        std::string line;
        // eat all leading comments
		while (is.good() && !is.eof() && fetch_line(is, line)) {
            // trim the line from leading and trailing space
            boost::trim(line);
            // skip empty lines
            if (line.empty())
                continue;
            // skip comments
            if (line[0] == '#')
               continue;

            // tokenize the line
            tokenizer tokens(line, sep);

            // access the tokens through a vector
            std::vector<std::string> vec;
            vec.assign(tokens.begin(),tokens.end());
            // get number of tohens
            int nrOfTokens = vec.size();
            /* Substract two from torken count and divide by 10 to get
             * the number of elements in the current span ...
             */
            int nrOfElements = (nrOfTokens - 2) / 10;
            /* ...and make sure the result (saved in an integer)
             * is the same as doing this backwards calculated witk
             * the tokens size. If not, complain in getting a wrong
             * size. This makes sure we have 12, 22, 32, ... tokens. The
             * patter is: having the map position as x and y followed by
             * 10 numbers per height span element.
             */
            if (nrOfTokens != 2 + 10 * nrOfElements) {
                std::stringstream err;
                int elements = (nrOfTokens - 2) % 10;
                err <<  "Error while loading version 2 file: Incomplete element data (got: "
                    << elements << " expected: 10)";
                Log(EWarn, err.str().c_str());
                continue;
            }

            int x;
            try {
                x = boost::lexical_cast<int>(vec[0]);
            } catch(boost::bad_lexical_cast &) {
                Log(EWarn, "Could not parse x cell position");
                continue;
            }

            int y;
            try {
                y = boost::lexical_cast<int>(vec[1]);
            } catch(boost::bad_lexical_cast &) {
                Log(EWarn, "Could not parse y cell position");
                continue;
            }

            // gat current heigh span
            HeightSpanType& hs = at(x,y);
            // iterate over its elements
            for (int i=0; i<nrOfElements; ++i) {
                float h1, h2;
			    int dists[4], indices[4];
                try {
                    int offset = 10 * i;
                    h1 = boost::lexical_cast<float>(vec[offset + 2]);
                    h2 = boost::lexical_cast<float>(vec[offset + 3]);

                    /* Read in the distances and indices to/of the adjacent spans.
                     * Make sure we get all four directions complete. If this is not
                     * possible, skip this element.
                     */
                    int j;
                    for (j = 0; j < 4; ++j) {
                        dists[j] = boost::lexical_cast<int>(vec[offset + 4 + 2*j]);
                        indices[j] = boost::lexical_cast<int>(vec[offset + 5 + 2*j]);
                    }
                    /* Everything seems alright, create the new element and put
                     * it into the collection of height span elements of the current
                     * cell.
                     */
                    hs.push_back( HeightSample(h1, h2) );
                    // iterate over neighbor directions
                    for (j = 0; j < 4; ++j) {
                        // index mapping from data format to internal structure
                        static int idx_trans[4] = { 2, 3, 0, 1 };
                        // get current internal structure index
                        int k = idx_trans[j];
                        // check if distance in current diretion (file) is zero
                        if (dists[j] == 0) {
                            /* Yes, so remember one as distance in that direction
                             * of internal structure. And remember -1 as slab index
                             * in that direction.
                             */
                            hs.back().neighbor_distances[k]    = 1;
                            hs.back().neighbor_slab_indices[k] = -1;
                        } else {
                            /* No, flag current direction to state existence of
                             * neigbor. Remember actual distance in that direction
                             * and the index of the adjacent span as well.
                             */
                            hs.back().set_flag(k);
                            hs.back().neighbor_distances[k]    = dists[j];
                            hs.back().neighbor_slab_indices[k] = indices[j];
                        }
                    }
                } catch(boost::bad_lexical_cast &) {
                    Log(EWarn, "Could not parse element data");
                    continue;
                }
            }
        }

        // do a sanity check
        check_connectivity();
    }

    void readDataV3(fs::ifstream& is) {
        separator sep("\t ");
        std::string line;
        // eat all leading comments
		while (is.good() && !is.eof() && fetch_line(is, line)) {
            // trim the line from leading and trailing space
            boost::trim(line);
            // skip empty lines
            if (line.empty())
                continue;
            // skip comments
            if (line[0] == '#')
               continue;

            // tokenize the line
            tokenizer tokens(line, sep);

            // access the tokens through a vector
            std::vector<std::string> vec;
            vec.assign(tokens.begin(),tokens.end());
            // get number of tohens
            int nrOfTokens = vec.size();
            /* Substract two from torken count and divide by 13 to get
             * the number of elements in the current span ...
             */
            int nrOfElements = (nrOfTokens - 2) / 13;
            /* ...and make sure the result (saved in an integer)
             * is the same as doing this backwards calculated witk
             * the tokens size. If not, complain in getting a wrong
             * size. This makes sure we have 15, 28, 41, ... tokens. The
             * patter is: having the map position as x and y followed by
             * 13 numbers per height span element.
             */
            if (nrOfTokens != 2 + 13 * nrOfElements) {
                std::stringstream err;
                int elements = (nrOfTokens - 2) % 13;
                err <<  "Error loading version 3 file: Incomplete element data (got: "
                    << elements << " expected: 13)";
                Log(EWarn, err.str().c_str());
                continue;
            }

            int x;
            try {
                x = boost::lexical_cast<int>(vec[0]);
            } catch(boost::bad_lexical_cast &) {
                Log(EWarn, "Could not parse x cell position");
                continue;
            }

            int y;
            try {
                y = boost::lexical_cast<int>(vec[1]);
            } catch(boost::bad_lexical_cast &) {
                Log(EWarn, "Could not parse y cell position");
                continue;
            }

            // gat current heigh span
            HeightSpanType& hs = at(x,y);
            // iterate over its elements
            for (int i=0; i<nrOfElements; ++i) {
                double h1, h2, dx, dy;
			    int type, dists[4], indices[4];
                try {
                    int offset = 13 * i;
                    h1 = boost::lexical_cast<double>(vec[offset + 2]);
                    type = boost::lexical_cast<int>(vec[offset + 3]);
                    h2 = boost::lexical_cast<double>(vec[offset + 4]);
                    dx = boost::lexical_cast<double>(vec[offset + 5]);
                    dy = boost::lexical_cast<double>(vec[offset + 6]);

                    /* Read in the distances and indices to/of the adjacent spans.
                     * Make sure we get all four directions complete. If this is not
                     * possible, skip this element.
                     */
                    int j;
                    for (j = 0; j < 4; ++j) {
                        dists[j] = boost::lexical_cast<int>(vec[offset + 7 + 2*j]);
                        indices[j] = boost::lexical_cast<int>(vec[offset + 8 + 2*j]);
                    }


                    /* Everything seems alright, create the new element and put
                     * it into the collection of height span elements of the current
                     * cell.
                     */
                    hs.push_back(HeightSample((float)h1,type,(float)h2,(float)dx,(float)dy));
                    // iterate over neighbor directions
                    for (j = 0; j < 4; ++j) {
                        // index mapping from data format to internal structure
                        static int idx_trans[4] = { 2, 3, 0, 1 };
                        // get current internal structure index
                        int k = idx_trans[j];
                        // check if distance in current diretion (file) is zero
                        if (dists[j] == 0) {
                            /* Yes, so remember one as distance in that direction
                             * of internal structure. And remember -1 as slab index
                             * in that direction.
                             */
                            hs.back().neighbor_distances[k]    = 1;
                            hs.back().neighbor_slab_indices[k] = -1;
                        } else {
                            /* No, flag current direction to state existence of
                             * neigbor. Remember actual distance in that direction
                             * and the index of the adjacent span as well.
                             */
                            hs.back().set_flag(k);
                            hs.back().neighbor_distances[k]    = dists[j];
                            hs.back().neighbor_slab_indices[k] = indices[j];
                        }
                    }
                } catch(boost::bad_lexical_cast &) {
                    Log(EWarn, "Could not parse element data");
                    continue;
                }
            }
        }

        // do a sanity check
        check_connectivity();
    }

    HeightSpanMap(Stream *stream, InstanceManager *manager) : Shape(stream, manager) {
		m_aabb = AABB(stream);
		m_name = stream->readString();
		unsigned int meshCount = stream->readUInt();
		m_meshes.resize(meshCount);

		for (unsigned int i=0; i<meshCount; ++i) {
			m_meshes[i] = static_cast<TriMesh *>(manager->getInstance(stream));
			m_meshes[i]->incRef();
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);

		m_aabb.serialize(stream);
		stream->writeString(m_name);
		stream->writeUInt((unsigned int) m_meshes.size());
		for (size_t i=0; i<m_meshes.size(); ++i)
			manager->serialize(stream, m_meshes[i]);
	}

	virtual ~HeightSpanMap() {
		for (size_t i=0; i<m_meshes.size(); ++i)
			m_meshes[i]->decRef();
		for (std::map<std::string, BSDF *>::iterator it = m_materials.begin();
			it != m_materials.end(); ++it) {
			(*it).second->decRef();
		}
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();
		if (cClass->derivesFrom(BSDF::m_theClass)) {
			m_bsdf = static_cast<BSDF *>(child);
			for (size_t i=0; i<m_meshes.size(); ++i) 
				m_meshes[i]->addChild(name, child);
			Assert(m_meshes.size() > 0);
			m_bsdf->setParent(NULL);
		} else if (cClass->derivesFrom(Luminaire::m_theClass)) {
			Assert(m_luminaire == NULL && m_meshes.size() == 1);
			m_luminaire = static_cast<Luminaire *>(child);
			for (size_t i=0; i<m_meshes.size(); ++i) {
				child->setParent(m_meshes[i]);
				m_meshes[i]->addChild(name, child);
			}
		} else if (cClass->derivesFrom(Subsurface::m_theClass)) {
			Assert(m_subsurface == NULL);
			m_subsurface = static_cast<Subsurface *>(child);
			for (size_t i=0; i<m_meshes.size(); ++i) { 
				child->setParent(m_meshes[i]);
				m_meshes[i]->addChild(name, child);
			}
		} else {
			Shape::addChild(name, child);
		}
	}

	bool isCompound() const {
		return true;
	}

	Shape *getElement(int index) {
		if (index >= (int) m_meshes.size())
			return NULL;
		Shape *shape = m_meshes[index];
		BSDF *bsdf = shape->getBSDF();
		Luminaire *luminaire = shape->getLuminaire();
		Subsurface *subsurface = shape->getSubsurface();
		if (bsdf)
			bsdf->setParent(shape);
		if (luminaire)
			luminaire->setParent(shape);
		if (subsurface)
			subsurface->setParent(shape);
		return shape;
	}

	std::string getName() const {
		return m_name;
	}

	AABB getAABB() const {
		return m_aabb;
	}

	Float getSurfaceArea() const {
		Float sa = 0;
		for (size_t i=0; i<m_meshes.size(); ++i)
			sa += m_meshes[i]->getSurfaceArea();
		return sa;
	}

	struct Vertex {
		Point p;
		Normal n;
		Point2 uv;
	};

	/// For using vertices as keys in an associative structure
	struct vertex_key_order : public 
		std::binary_function<Vertex, Vertex, bool> {
	public:
		bool operator()(const Vertex &v1, const Vertex &v2) const {
			if (v1.p.x < v2.p.x) return true;
			else if (v1.p.x > v2.p.x) return false;
			if (v1.p.y < v2.p.y) return true;
			else if (v1.p.y > v2.p.y) return false;
			if (v1.p.z < v2.p.z) return true;
			else if (v1.p.z > v2.p.z) return false;
			if (v1.n.x < v2.n.x) return true;
			else if (v1.n.x > v2.n.x) return false;
			if (v1.n.y < v2.n.y) return true;
			else if (v1.n.y > v2.n.y) return false;
			if (v1.n.z < v2.n.z) return true;
			else if (v1.n.z > v2.n.z) return false;
			if (v1.uv.x < v2.uv.x) return true;
			else if (v1.uv.x > v2.uv.x) return false;
			if (v1.uv.y < v2.uv.y) return true;
			else if (v1.uv.y > v2.uv.y) return false;
			return false;
		}
	};

    /**
     * Initialize the mitsube usable geometry.
     */
    void generateGeometry(const std::string &name,
			bool hasNormals, bool hasTexcoords,
			BSDF *currentMaterial,
            const Transform &objectToWorld) {
		Log(EInfo, "Loading geometry \"%s\"", name.c_str());

		std::map<Vertex, int, vertex_key_order> vertexMap;
		std::vector<Vertex> vertexBuffer;
		size_t numMerged = 0;
		Vector translate(0.0f);
		Float scale = 0.0f;
        // index mapping for inverse directions
        static int idx_trans[4] = { 2, 3, 0, 1 };

		if (m_recenter) {
			AABB aabb;
            for (HeightSpanIterator hsi=begin(); !at_end(hsi); next(hsi)) {
				aabb.expandBy(getSurfaceSample(hsi));
            }
			scale = 2/aabb.getExtents()[aabb.getLargestAxis()];
			translate = -Vector(aabb.getCenter());
		}

        std::vector<Triangle> triangles;

        // make sure no previous user flag is set anymore
        clearUserFlag();
		/* Collapse the mesh into a more usable form */
        // iterate over all height spans/elements
        for (HeightSpanIterator hsi=begin(); !at_end(hsi); next(hsi)) {
            // get a reference to the current cell
            HeightSample& hs = refSample(hsi);
            // remember this vertex
            //vertices.push_back( Point() );
            // iterate over all four directions
            for (int l1=0; l1<4; ++l1) {
                // calculate the second point by rotating l1
                int l2 = (l1+1)&3;
                /* Check  if there actually is a tringle and if this
                 * triangle has not been touched already.
                 */
                if (hs.get_flag(l1) && hs.get_flag(l2) && !hs.get_user_flag(l1)) {
                    // alright, use this triangle and note this with a user flag
                    hs.set_user_flag(l1);
                    /* Get two new iterators on the points ending the
                     * current directions.
                     */
                    HeightSpanIterator hsi1 = hsi;
                    move(hsi1, l1);
                    HeightSpanIterator hsi2 = hsi;
                    move(hsi2, l2);
                    // do an inverse reference test
                    const HeightSample& hs1 = getSample(hsi1);
                    const HeightSample& hs2 = getSample(hsi2);
                    if (!hs1.get_flag( idx_trans[l1] ) || !hs2.get_flag( idx_trans[l2])
                            || (hs1.get_nbr_slab_idx( idx_trans[l1] ) != hsi.k)
                            || (hs2.get_nbr_slab_idx( idx_trans[l2] ) != hsi.k)) {
                        Log(EWarn, "Ignoring triangle because of missing inverse connection information");
                        continue;
                    }
                    // calculate the actual positions (still in object space)
                    std::vector<Point> pts(3);
                    pts[0] = getSurfaceSample(hsi);
                    pts[1] = getSurfaceSample(hsi1);
                    pts[2] = getSurfaceSample(hsi2);
                    Normal n = cross(pts[2] - pts[0], pts[1] - pts[0]);

                    // create triangle
                    Triangle tri;

                    for (unsigned int j=0; j<3; ++j) {
                        int key;        
                        Vertex vertex;
 
                        if (m_recenter)
                            vertex.p = objectToWorld( m_fileToObject( (pts[j] + translate) * scale ) );
                        else
                            vertex.p = objectToWorld( m_fileToObject( pts[j] ) );

                        vertex.n = Normal(0.0f);
					    vertex.uv = Point2(0.0f);

                        if (vertexMap.find(vertex) != vertexMap.end()) {
                            key = vertexMap[vertex];
                            numMerged++;
                        } else {
                            key = (int) vertexBuffer.size();
                            vertexMap[vertex] = (int) key;
                            vertexBuffer.push_back(vertex);
                        }

                        tri.idx[j] = key;
                    }
                    triangles.push_back(tri);
                    
                }
            }
        }		
		
        ref<TriMesh> mesh = new TriMesh(name,
			triangles.size(), vertexBuffer.size(),
			hasNormals, hasTexcoords, false,
			m_flipNormals, m_faceNormals);

		std::copy(triangles.begin(), triangles.end(), mesh->getTriangles());

		Point    *target_positions = mesh->getVertexPositions();
		Normal   *target_normals   = mesh->getVertexNormals();
		Point2   *target_texcoords = mesh->getVertexTexcoords();

		for (size_t i=0; i<vertexBuffer.size(); i++) {
			*target_positions++ = vertexBuffer[i].p;
			if (hasNormals)
				*target_normals++ = vertexBuffer[i].n;
			if (hasTexcoords)
				*target_texcoords++ = vertexBuffer[i].uv;
		}

		mesh->incRef();
		if (currentMaterial)
			mesh->addChild("", currentMaterial);
		m_meshes.push_back(mesh);
		Log(EInfo, "%s: Loaded " SIZE_T_FMT " triangles, " SIZE_T_FMT 
			" vertices (merged " SIZE_T_FMT " vertices).", name.c_str(),
			triangles.size(), vertexBuffer.size(), numMerged);
		mesh->configure();
    }
	
	void configure() {
		Shape::configure();

		m_aabb.reset();
		for (size_t i=0; i<m_meshes.size(); ++i) {
			m_meshes[i]->configure();
			m_aabb.expandBy(m_meshes[i]->getAABB());
		}
	}

    /**
     * Remove any user flags that might have been set.
     */
    void clearUserFlag()
    {
        for (int i=0; i<m_width; ++i)
            for (int j=0; j<m_height; ++j) {
                HeightSpanType& hst = at(i,j);
                if (hst.empty())
                    continue;
                for (int k=0; k < (int)hst.size(); ++k)
                    hst[k].user_flags = 0;
            }
    }

    /**
     * Resize map and set size members
     */
	void resize(int w, int h)
	{
		m_width = w;
		m_height = h;
		m_heightSpans.resize(w * h);
	}

    /**
     * Returns height span at position (i,j).
     */
	HeightSpanType& at(int i, int j) {
        return m_heightSpans[j * m_width + i];
    }

    /**
     * Returns const height span at position (i,y).
     */
    const HeightSpanType& at(int i, int j) const {
        return m_heightSpans[j * m_width + i];
    }

    /**
     * Get map index diffence between two adjacent cells in x direction.
     */
	static int getDeltaX(int i) {
		static int dx[4] = { 0, -1, 0, 1 };
		return dx[i];
	}

    /**
     * Get map index diffence between two adjacent cells in y direction.
     */
	static int getDeltaY(int i) {
		static int dy[4] = { -1, 0, 1, 0 };
		return dy[i];
	}

    /**
     * Get direction opposite to direction i.
     */
	static int getNeighborIndex(int i) {
		static int ni[4] = { 2, 3, 0, 1 };
		return ni[i];
	}

    /**
     * Checks if the index pair is inside the map.
     */
	bool isInside(int x, int y) const {
		return x >= 0 && y >= 0 && x < m_width && y < m_height;
	}

    /**
     * Get a the current element from the current iterator position (i.e. current height span).
     */
	const HeightSample& getSample(const HeightSpanIterator& hsi) const {
         return at(hsi.i, hsi.j)[hsi.k];
    }

    /**
     * Get current element of the passed iterator.
     */
	HeightSample& refSample(const HeightSpanIterator& hsi) {
        return at(hsi.i, hsi.j)[hsi.k];
    }

    /**
     * Get the 3D position of the top of the iterators current element.
     */
    Point getSurfaceSample(const HeightSpanIterator& hsi) const
    {
        const HeightSample& hs = getSample(hsi);
        return Point( hsi.x + hsi.dx * hs.dx, hsi.y + hsi.dy * hs.dy, hs.get_total_height());
    }

    /**
     * Moves the iterator hsi one step in direction l. Returns true
     * on success and false otherwise (e.g. on border).
     */
	bool move(HeightSpanIterator& hsi, int l) const {
        // get current element of iterator
        const HeightSample& hs = getSample(hsi);
        // return unsuccessful if there is no neighbor in in direction l
        if (!hs.get_flag(l))
            return false;
        // get distance of current element to neighbor in direction l
        int dist = hs.get_distance(l);
        // calculate index of neighbor in direction l
        int n_i = hsi.i + dist * getDeltaX(l);
        int n_j = hsi.j + dist * getDeltaY(l);
        // make sure the neighbor is inside the map
        if (!isInside(n_i,n_j))
            return false;
        // move the iterator to the neighbor position
        hsi.i  = n_i;
        hsi.j  = n_j;
        hsi.x += hsi.dx * dist * getDeltaX(l);
        hsi.y += hsi.dy * dist * getDeltaY(l);
        hsi.k  = hs.get_nbr_slab_idx(l);
        return true;
    }


    /**
     * Move the passed iterato to the next element, potentially in
     * another cell. Will return true if successful, false otherwise
     * (i.e. no more elements are on this map after the current one).
     */
	bool next(HeightSpanIterator& hsi) const {
        // get current height span of iterator
        const HeightSpanType& hst = at(hsi.i, hsi.j);
        // check if there is a next element in current height span
        if (hsi.k + 1 < (int)hst.size()) {
            // yes, so increment offset within span and return successful
            ++hsi.k;
            return true;
        }
        // no, so start with zero offset and look for next span
        hsi.k = 0;
        while (true) {
            // check if there is next span in height direction
            if (hsi.j+1 == m_height) {
                // no, reset height position information in iterator
                hsi.j = 0;
                hsi.y = 0; //0.5f*hsi.dy;
                // increment width position information by span distance
                hsi.x += hsi.dx;
                /* Increment and check if there is a next span in width
                 * direction.
                 */
                if (++hsi.i == m_width)
                    // no, return unsuccessful
                    return false;
            } else {
                // yes, go to next span in height direction
                ++hsi.j;
                hsi.y += hsi.dy;
            }
            // check if the current (just yet arrived at) span is not empy
            if (!at(hsi.i, hsi.j).empty())
                // it is not empty, return unsuccessful
                return true;
            // go to next span if not yet returned
        }
    }

    /**
     * Get an iterator positioned at the first element of the first
     * span of the map.
     */
    HeightSpanIterator begin() const {
        int i,j;
        bool found = false;
        /* Walk the map, until the first heigt span with at least
         * one element is found.
         */
        for (i=0; !found && i<m_width; ++i)
            for (j=0; !found && j<m_height; ++j) {
                const HeightSpanType& hst = at(i,j);
                // are some elements in this cell?
                if (!hst.empty())
                    // yes, cool, return this cell/height span
                    return HeightSpanIterator(1.0f/(m_width-1), 1.0f/(m_height-1), i, j);
            }
        // no height spans fourd and hence no snow found, complain!
        Log(EError, "Oops, no height sample found in map");
        // return dummy iterator
        return HeightSpanIterator(1.0f/(m_width-1), 1.0f/(m_height-1), -1, -1);
    }

    /**
     * Check if the iterator is after the last element, i.e. the end.
     */
    bool at_end(const HeightSpanIterator& hsi) const
    {
        // we are at the end if the width index is out of range
        return hsi.i == m_width;
    }

    /**
     *  Does a sanity check for the interconnection of elemennts.
     */
    void check_connectivity() {
        /* start with first element of first height span and walk
         * elementwise through the map.
         */
        for (HeightSpanIterator hsi=begin(); !at_end(hsi); next(hsi)) {
            // get the current element of the current iterator position
            const HeightSample& hs = getSample(hsi);
            // iterate over all four neigbor directions
            for (int l=0; l<4; ++l) {
                // check if there is a neighbor in the current direction
                if (hs.get_flag(l)) {
                    // yes, check for its index
                    if (hs.get_nbr_slab_idx(l) == -1) {
                        // inconsistency detected, no valid index found
                        Log(EWarn, "Inconsistend nbr slab index and flag");
                        continue;
                    }
                    // A distance below one is not allowed to neighbors, check for it
                    if (hs.get_distance(l) < 1) {
                        Log(EWarn, "Invalid neighbor distance");
                        continue;
                    }
                    // get indices of current neighbor
                    int n_i = hsi.i + getDeltaX(l) * hs.get_distance(l);
                    int n_j = hsi.j + getDeltaY(l) * hs.get_distance(l);
                    // check if it is inside map index bounds.
                    if (!isInside(n_i, n_j)) {
                        // it is not, complain!
                        Log(EWarn, "reference to neighbor outside of valid range");
                        continue;
                    }
                    // get actual neigboring element
                    HeightSpanType& n_hst = at(n_i,n_j);
                    /* Check if our slab reference is larger than number of
                     * slabs in referenced cell.
                     */
                    if (hs.get_nbr_slab_idx(l) >= (int)n_hst.size()) {
                        // it is larger, complain!
                        std::stringstream err;
                        err << "(" << hsi.i << "," << hsi.j << "): reference to (dir: " << l << ") neighbor ("
                            << n_i << "," << n_j << ") slab " << hs.get_nbr_slab_idx(l)
                            << " outside of valid range [0," << n_hst.size() << "[";
                        Log(EWarn, err.str().c_str());
                        continue;
                    }
                }
            }
        }
    }

	MTS_DECLARE_CLASS()
private:
    /* members */
	std::vector<TriMesh *> m_meshes;
	std::map<std::string, BSDF *> m_materials;
	bool m_flipNormals, m_faceNormals, m_recenter;
	std::string m_name;
	AABB m_aabb;
    // dimensions of the map
    int m_width, m_height;
    // the actual map, saved in a linear collection.
	std::vector<HeightSpanType> m_heightSpans;
    /* the transformation used to convert file data
     * to actuel 3d space data.
     */
     Transform m_fileToObject;
};

MTS_IMPLEMENT_CLASS_S(HeightSpanMap, false, Shape)
MTS_EXPORT_PLUGIN(HeightSpanMap, "Height Span Map loader");
MTS_NAMESPACE_END
