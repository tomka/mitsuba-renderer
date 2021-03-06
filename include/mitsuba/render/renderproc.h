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

#if !defined(__RENDERPROC_H)
#define __RENDERPROC_H

#include <mitsuba/render/scene.h>
#include <mitsuba/render/imageproc.h>
#include <mitsuba/render/renderqueue.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Parallel process for rendering with sampling-based integrators.
 *
 * Splits an image into independent rectangular pixel regions, which are
 * then rendered in parallel.
 *
 * \sa SampleIntegrator
 */
class MTS_EXPORT_RENDER BlockedRenderProcess : public BlockedImageProcess {
public:
	BlockedRenderProcess(const RenderJob *parent, RenderQueue *queue, 
		int blockSize);

	// ======================================================================
	//! @{ \name Implementation of the ParallelProcess interface
	// ======================================================================

	ref<WorkProcessor> createWorkProcessor() const;
	void processResult(const WorkResult *result, bool cancelled);
	void bindResource(const std::string &name, int id);
	EStatus generateWork(WorkUnit *unit, int worker);
	
	//! @}
	// ======================================================================

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~BlockedRenderProcess();
protected:
	ref<RenderQueue> m_queue;
	ref<Scene> m_scene;
	ref<Film> m_film;
	const RenderJob *m_parent;
	int m_resultCount;
	ref<Mutex> m_resultMutex;
	ProgressReporter *m_progress;
	int m_borderSize;
};

MTS_NAMESPACE_END

#endif /* __RENDERPROC_H */
