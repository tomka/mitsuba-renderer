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

#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/mipmap.h>

MTS_NAMESPACE_BEGIN

/**
 * Simple linear (i.e. not gamma corrected) bitmap texture
 * using the EXR file format
 */
class DiffusionTexture : public Texture2D {
public:
	DiffusionTexture(const Properties &props) : Texture2D(props) {
		std::string filterType = props.getString("filterType", "ewa");
		std::string wrapMode = props.getString("wrapMode", "repeat");

		if (filterType == "ewa")
			m_anisotropic = true;
		else if (filterType == "isotropic")
			m_anisotropic = false;

		if (wrapMode == "repeat")
			m_wrapMode = MIPMap::ERepeat;
		else if (wrapMode == "clamp")
			m_wrapMode = MIPMap::EClamp;
		else if (wrapMode == "black")
			m_wrapMode = MIPMap::EBlack;
		else if (wrapMode == "white")
			m_wrapMode = MIPMap::EWhite;
		else
			Log(EError, "Unknown wrap mode '%s' -- must be "
				"'repeat', 'clamp', 'black', or 'white'!", filterType.c_str());

		m_maxAnisotropy = props.getFloat("maxAnisotropy", 8);

        bool loadDataManually = props.getBoolean("loadDataManually", false);

        if (!loadDataManually) {
            m_filename = Thread::getThread()->getFileResolver()->resolve(
                props.getString("filename"));
            Log(EInfo, "Loading texture \"%s\"", m_filename.leaf().c_str());

            ref<FileStream> fs = new FileStream(m_filename, FileStream::EReadOnly);
        
            ref<Bitmap> bitmap = new Bitmap(Bitmap::EEXR, fs);
            initializeFrom(bitmap, false);
        }
	}

	DiffusionTexture(Stream *stream, InstanceManager *manager) 
	 : Texture2D(stream, manager) {
		m_filename = stream->readString();
		Log(EInfo, "Unserializing texture \"%s\"", m_filename.leaf().c_str());
		stream->writeBool(m_anisotropic);
		stream->writeUInt(m_wrapMode);
		stream->writeFloat(m_maxAnisotropy);
		int size = stream->readInt();
		ref<MemoryStream> mStream = new MemoryStream(size);
		stream->copyTo(mStream, size);
		mStream->setPos(0);
		ref<Bitmap> bitmap = new Bitmap(Bitmap::EEXR, mStream);
        initializeFrom(bitmap, false);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture2D::serialize(stream, manager);
		stream->writeString(m_filename.file_string());
		stream->writeBool(m_anisotropic);
stream->writeUInt(m_wrapMode);
		stream->writeFloat(m_maxAnisotropy);
		ref<Stream> is = new FileStream(m_filename, FileStream::EReadOnly);
		stream->writeInt(is->getSize());
		is->copyTo(stream);
	}

	void initializeFrom(Bitmap *bitmap, bool copy) {
        if (copy) {
    		ref<Bitmap> own = new Bitmap(bitmap->getWidth(), bitmap->getHeight(), 128);

            float *data = bitmap->getFloatData();
            float *flData = own->getFloatData();

            for (int y=0; y<bitmap->getHeight(); ++y) {
                for (int x=0; x<bitmap->getWidth(); ++x) {
                    *flData++ = *data++;
                    *flData++ = *data++;
                    *flData++ = *data++;
                    *flData++ = *data++;
                }
            }
            bitmap = own;
        }

		// m_mipmap = MIPMap::fromBitmap(bitmap);
		m_mipmap = MIPMap::fromBitmap(bitmap, !m_anisotropic,
				m_wrapMode, m_maxAnisotropy);
		m_average = m_mipmap->triangle(m_mipmap->getLevels()-1, 0, 0);
		m_maximum = m_mipmap->getMaximum();
    }

	Spectrum getValue(const Point2 &uv) const {
		return m_mipmap->triangle(0, uv.x, uv.y);
	}

	Spectrum getValue(const Point2 &uv, Float dudx, 
			Float dudy, Float dvdx, Float dvdy) const {
		return m_mipmap->getValue(uv.x, uv.y, dudx, dudy, dvdx, dvdy);
	}

	Spectrum getMaximum() const {
		return m_maximum;
	}

	Spectrum getAverage() const {
		return m_average;
	}

	bool usesRayDifferentials() const {
		return true;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "DiffusionTexture[filename=\"" << m_filename.file_string() << "\"]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	ref<MIPMap> m_mipmap;
	fs::path m_filename;
	Spectrum m_average, m_maximum;
	bool m_anisotropic;
	MIPMap::EWrapMode m_wrapMode;
	Float m_maxAnisotropy;
};

MTS_IMPLEMENT_CLASS_S(DiffusionTexture, false, Texture2D)
MTS_EXPORT_PLUGIN(DiffusionTexture, "HDR texture (EXR)");
MTS_NAMESPACE_END
