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

#include <mitsuba/core/timer.h>
#include <mitsuba/render/renderjob.h>

MTS_NAMESPACE_BEGIN

RenderQueue::RenderQueue(EExecutionStrategy execStrategy) {
    m_managingStrategy = execStrategy; 
	m_mutex = new Mutex();
	m_joinMutex = new Mutex();
	m_cond = new ConditionVariable(m_mutex);
	m_timer = new Timer();
}
	
RenderQueue::~RenderQueue() {
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->decRef();
}

void RenderQueue::addJob(RenderJob *job) {
	m_mutex->lock();
	m_jobs[job] = JobRecord(m_timer->getMilliseconds());
	job->incRef();
	m_mutex->unlock();
}

void RenderQueue::registerListener(RenderListener *listener) {
	listener->incRef();
	m_mutex->lock();
	m_listeners.push_back(listener);
	m_mutex->unlock();
}

void RenderQueue::unregisterListener(RenderListener *listener) {
	m_mutex->lock();
	m_listeners.erase(std::remove(m_listeners.begin(), m_listeners.end(), listener));
	m_mutex->unlock();
	listener->decRef();
}
	
void RenderQueue::flush() {
	m_mutex->lock();
	std::map<RenderJob *, JobRecord>::iterator it = m_jobs.begin();
	for (; it != m_jobs.end(); ++it) {
		(*it).first->flush();
	}
	m_mutex->unlock();
}

void RenderQueue::managedExecution(RenderJob *thr) {
	m_mutex->lock();
	std::map<RenderJob *, JobRecord>::iterator it = m_jobs.find(thr);
	if (it == m_jobs.end()) {
		Log(EError, "RenderQueue::managedExecution() - job not found!");
		m_mutex->unlock();
	}
	JobRecord &rec = (*it).second;
    if (rec.delayed)
		Log(EWarn, "RenderQueue::managedExecution() - job already set up for delayed execution!");

    // enqueue if this is not the first job
    bool queueJob = (m_jobs.size() > 1) && (m_managingStrategy == ESerial);
    if (queueJob) {
        m_waitingJobs.push(thr);
	    Log(EDebug, "Managing execution of new job: queued");
    } else {
        thr->start();
	    Log(EDebug, "Managing execution of new job: run directly");
    }
    m_mutex->unlock();
}

void RenderQueue::setManagedExecutionStrategy(EExecutionStrategy es) {
	m_mutex->lock();
	if (m_jobs.size() > 0) {
		Log(EWarn, "Job queue not empy - won't change managment strategy!");
        m_mutex->unlock();
        return;
	}

    m_managingStrategy = es;
	Log(EDebug, "Set new managed execution strategy");
    m_mutex->unlock();
}

void RenderQueue::removeJob(RenderJob *job, bool cancelled) {
	m_mutex->lock();
	std::map<RenderJob *, JobRecord>::iterator it = m_jobs.find(job);
	if (it == m_jobs.end()) {
		Log(EError, "RenderQueue::removeRenderJob() - job not found!");
		m_mutex->unlock();
	}
	JobRecord &rec = (*it).second;
	unsigned int ms = m_timer->getMilliseconds() - rec.startTime;
    if (rec.delayed)
	    Log(EInfo, "Render time: %s", timeString(ms/1000.0f, true).c_str());
    else
    	Log(EInfo, "Render time: %s after waiting %s", timeString(ms/1000.0f, true).c_str(),
            timeString(rec.waitTime/1000.0f, true).c_str());
	m_jobs.erase(job);
	m_cond->broadcast();
	m_joinMutex->lock();
	m_joinList.push_back(job);
	m_joinMutex->unlock();
	signalFinishJob(job, cancelled);

    // execute a potentially delayed job
    if (m_waitingJobs.size() > 0) {
        RenderJob *waitingJob = m_waitingJobs.front();
        m_waitingJobs.pop();

        it = m_jobs.find(waitingJob);
        if (it == m_jobs.end()) {
            Log(EError, "RenderQueue::removeRenderJob() - queued job not found!");
            m_mutex->unlock();
        }
	    JobRecord waitRec = (*it).second;
        // set start time
        unsigned int origStartTime = waitRec.startTime;
        waitRec.startTime = m_timer->getMilliseconds();
        waitRec.waitTime = m_timer->getMilliseconds() - origStartTime;
        m_jobs[waitingJob] = waitRec;
        // finally, start the waiting job
        waitingJob->start();
    }
	m_mutex->unlock();
}
	
void RenderQueue::waitLeft(size_t njobs) const {
	m_mutex->lock();
	while (m_jobs.size() > njobs) 
		m_cond->wait();
	m_mutex->unlock();
	join();
}

void RenderQueue::join() const {
	m_joinMutex->lock();
	/* Wait for the proper termination of all stopping threads */
	for (size_t i=0; i<m_joinList.size(); ++i) {
		RenderJob *job = m_joinList[i];
		job->join();
		job->decRef();
	}
	m_joinList.clear();
	m_joinMutex->unlock();
}

void RenderQueue::signalWorkBegin(const RenderJob *job, const RectangularWorkUnit *wu, int worker) {
	m_mutex->lock();
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->workBeginEvent(job, wu, worker);
	m_mutex->unlock();
}

void RenderQueue::signalWorkEnd(const RenderJob *job, const ImageBlock *wr) {
	m_mutex->lock();
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->workEndEvent(job, wr);
	m_mutex->unlock();
}

void RenderQueue::signalFinishJob(const RenderJob *job, bool cancelled) {
	m_mutex->lock();
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->finishJobEvent(job, cancelled);
	m_mutex->unlock();
}

void RenderQueue::signalRefresh(const RenderJob *job, const Bitmap *bitmap) {
	m_mutex->lock();
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->refreshEvent(job, bitmap);
	m_mutex->unlock();
}

MTS_IMPLEMENT_CLASS(RenderQueue, false, Object)
MTS_IMPLEMENT_CLASS(RenderListener, true, Object)
MTS_NAMESPACE_END
