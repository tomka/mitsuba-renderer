/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__SUBSURFACE_H)
#define __SUBSURFACE_H

#include <mitsuba/core/netobject.h>

MTS_NAMESPACE_BEGIN

/**
 * A storage place for material properties for the subsurface
 * integrators in use. It is implemented as a singlton object
 * and offers speed-ups when used e.g. for LUT storage.
 * ToDo: Add mutex.
 */
class MTS_EXPORT_RENDER SubsurfaceMaterialManager : public Object {

    /* private constructor */
    SubsurfaceMaterialManager( ) { }

public:
    typedef std::vector<Spectrum> LUTType;

    /**
     * Small data structure to store information about a look-up-table.
     */
    struct LUTRecord {
        Float resolution;
        ref<LUTType> lut;

        /* create invalid LUTRecord */
        LUTRecord()
                : resolution(-1), lut(NULL) { }

        LUTRecord(ref<LUTType> _lut, Float res)
                : resolution(res) {
            lut = _lut;
        }
    };

	MTS_DECLARE_CLASS()
public:
    bool hasLUT(const std::string &hash) const;

    void addLUT(const std::string &hash, const LUTRecord &lutRec);

    LUTRecord getLUT(const std::string &hash) const;

    std::string getMultipoleLUTHash(Float resolution, Float errorThreshold,
        const Spectrum &sigmaTr, const Spectrum &alphaPrime, int numExtraDipoles,
        const std::vector<Spectrum> &zrList, const std::vector<Spectrum> &zvList) const;

    std::string getDipoleLUTHash(Float resolution, Float errorThreshold,
        const Spectrum &sigmaTr, const Spectrum &alphaPrime, 
        const Spectrum &zr, const Spectrum &zv) const;

    static SubsurfaceMaterialManager* getInstance() {
        if (m_instance.get() == NULL) {
            m_instance = new SubsurfaceMaterialManager();
        }
        return m_instance;
    }

protected:
    typedef std::map<std::string, LUTRecord> StorageType;

    /* singlton instance */
    static ref<SubsurfaceMaterialManager> m_instance;

    /* lut map */
    StorageType m_lutRecords;

	/* Virtual destructor */
	virtual ~SubsurfaceMaterialManager();
};

/**
 * Abstract subsurface integrator -- can be attached to an arbitrary 
 * shape to compute exitant radiance due to internal scattering.
 */
class MTS_EXPORT_RENDER Subsurface : public NetworkedObject {
public:
	/**
	 * Possibly perform a pre-process task. The last three parameters are
	 * resource IDs of the associated scene, camera and sample generator,
	 * which have been made available to all local and remote workers.
	 */
	virtual bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID) = 0;

	/// Cancel any running pre-process tasks
	virtual void cancel();

	/// Return the list of shapes associated with this subsurface integrator
	inline const std::vector<Shape *> getShapes() const { return m_shapes; }

	/// Get the exitant radiance for a point on the surface
	virtual Spectrum Lo(const Scene *scene, Sampler *sampler,
		const Intersection &its, const Vector &d, int depth = 0) const = 0;

	/// Serialize this subsurface integrator to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	/// Set the parent node of the subsurface integrator
	void setParent(ConfigurableObject *parent);

	MTS_DECLARE_CLASS()
protected:
	/// Create a new subsurface scattering class
	Subsurface(const Properties &props);

	/// Unserialize a subsurface integrator from a binary data stream
	Subsurface(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~Subsurface();
protected:
	Spectrum m_sigmaS;
	Spectrum m_sigmaA;
	Spectrum m_sigmaT;
	Float m_eta, m_densityMultiplier;
	std::vector<Shape *> m_shapes;
};

MTS_NAMESPACE_END

#endif /* __SUBSURFACE_H */
