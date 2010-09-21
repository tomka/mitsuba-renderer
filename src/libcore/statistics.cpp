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

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/lock.h>

MTS_NAMESPACE_BEGIN

// -----------------------------------------------------------------------
//  Statistics collection
// -----------------------------------------------------------------------

bool ProgressReporter::m_enabled = true;

ProgressReporter::ProgressReporter(const std::string &title, long long total, const void *ptr)
 : m_title(title), m_total(total), m_value(0), m_percentage(-1), m_fillPos(0), m_ptr(ptr) {
	m_fillSize = (int) (PROGRESS_MSG_SIZE - title.length() - 3);
	SAssert(m_fillSize > 0);
	for (int i=0; i<m_fillSize; i++)
		m_string[i] = ' ';
	m_string[m_fillSize] = '\0';
	m_timer = new Timer();
	m_lastMs = 0;
}

void ProgressReporter::setEnabled(bool value) {
	m_enabled = value;
}

void ProgressReporter::reset() {
	for (int i=0; i<m_fillSize; i++)
		m_string[i] = ' ';
	m_timer->reset();
	m_lastMs = 0;
	m_value = 0;
	m_percentage = -1;
	m_fillPos = 0;
}

void ProgressReporter::update(long long value) {
	if (!m_enabled)
		return;
	value = std::min(std::max(value, (long long) 0), (long long) m_total);

	Float perc = (value * 100.0f) / m_total;
	unsigned int curMs = m_timer->getMilliseconds();
	m_value = value;

	if (value == m_total || (
			curMs - m_lastMs > 1000)) {
		m_percentage = (int) perc;
		int fillEnd = (int) ((value * m_fillSize) / m_total);
	
		Float time = curMs / 1000.0f;
		Float remaining = (time*m_total) / value - time;

		for (int i=m_fillPos; i<fillEnd; i++)
			m_string[i] = '+';
		m_fillPos = fillEnd;
		std::ostringstream oss;
		std::string eta = timeToString(remaining);
		oss << '\r' << m_title << ": [" << m_string << "] (";
		oss << timeToString(time) << ", ETA: " 
			<< eta << ")  \b\b";
		Thread::getThread()->getLogger()->logProgress(
			perc, m_title, oss.str(), eta, m_ptr);
		m_lastMs = curMs;
	}
}

StatsCounter::StatsCounter(const std::string &cat, const std::string &name, EStatsType type, uint64_t initial, uint64_t base)
 : m_category(cat), m_name(name), m_type(type) {
	m_value = (CacheLineCounter *) allocAligned(sizeof(CacheLineCounter) * NUM_COUNTERS);
	m_base = (CacheLineCounter *) allocAligned(sizeof(CacheLineCounter) * NUM_COUNTERS);
	memset(m_value, 0, sizeof(CacheLineCounter) * NUM_COUNTERS);
	memset(m_base, 0, sizeof(CacheLineCounter) * NUM_COUNTERS);
#if defined(WIN32) && !defined(WIN64)
	if (sizeof(LONG) != sizeof(uint32_t)) {
		cerr << "Internal error: sizeof(LONG) != sizeof(uint32_t)!" << endl;
		exit(-1);
	}
	m_value[0].value = (LONG) initial;
	m_base[0].value = (LONG) base;
#else
	m_value[0].value = initial;
	m_base[0].value = base;
#endif
	Statistics::getInstance()->registerCounter(this);
}

StatsCounter::~StatsCounter() {
	freeAligned(m_value);
	freeAligned(m_base);
}

bool StatsCounter::operator<(const StatsCounter &v) const {
	if (getCategory() == v.getCategory())
		return getName() < v.getName();
	return getCategory() < v.getCategory();
}

ref<Statistics> Statistics::m_instance = new Statistics();
	
void Statistics::staticInitialization() {
	SAssert(sizeof(CacheLineCounter) == 128);
}

void Statistics::staticShutdown() {
	m_instance = NULL;
}

Statistics::Statistics() {
	m_mutex = new Mutex();
}

void Statistics::registerCounter(const StatsCounter *ctr) {
	m_counters.push_back(ctr);
}

void Statistics::logPlugin(const std::string &name, const std::string &descr) {
	m_plugins.push_back(std::pair<std::string, std::string>(name, descr));
}

void Statistics::printStats() {
	SLog(EInfo, "Statistics: \n%s", getStats().c_str());
}

std::string Statistics::getStats() {
	std::ostringstream oss;
	m_mutex->lock();
	oss << "------------------------------------------------------------" << std::endl;

	oss << " * Loaded plugins :" << std::endl;

	if (m_plugins.size() == 0) {
		oss << "     none." << std::endl;
	} else {
		std::sort(m_plugins.begin(), m_plugins.end());
		for (unsigned int i=0; i<m_plugins.size(); i++) 
			oss << "    -  " << m_plugins[i].first << " [" << m_plugins[i].second
				<< "]" << endl;
	}

	std::string suffixesNumber[] = { "", " K", " M", " G", " T" };
	std::string suffixesByte[] = { " B", " KiB", " MiB", " GiB", " TiB" };
	std::string category = "";

	std::sort(m_counters.begin(), m_counters.end(), compareCategory());
	int statsEntries = 0;

	for (size_t i=0; i<m_counters.size(); ++i) {
		const StatsCounter *counter = m_counters[i];
		char temp[128];
		float value = (float) counter->getValue();
		float baseValue = (float) counter->getBase();
		int suffixIndex = 0;

		if (value == 0 && counter->getType() != EPercentage)
			continue;

		if (category != counter->getCategory()) {
			category = counter->getCategory();
			oss << std::endl <<  "  * " << category << " :" << std::endl;
		}

		switch (counter->getType()) {
			case ENumberValue:
				while (value > 1000.0f) {
					value /= 1000.0f;
					suffixIndex++;
				}
				if (value - std::floor(value) < 0.001f)
					snprintf(temp, sizeof(temp), "    -  %s : %.0f%s", counter->getName().c_str(), value, suffixesNumber[suffixIndex].c_str());
				else
					snprintf(temp, sizeof(temp), "    -  %s : %.3f%s", counter->getName().c_str(), value, suffixesNumber[suffixIndex].c_str());
				break;
			case EByteCount:
				while (value > 1024.0f) {
					value /= 1024.0f;
					suffixIndex++;
				}
				if (value - std::floor(value) < 0.001f)
					snprintf(temp, sizeof(temp), "    -  %s : %.0f%s", counter->getName().c_str(), value, suffixesByte[suffixIndex].c_str());
				else
					snprintf(temp, sizeof(temp), "    -  %s : %.3f%s", counter->getName().c_str(), value, suffixesByte[suffixIndex].c_str());
				break;
			case EPercentage: {
					Float value2 = value, value3 = baseValue;
					while (value2 > 1000.0f) {
						value2 /= 1000.0f;
						value3 /= 1000.0f;
						suffixIndex++;
					}
					snprintf(temp, sizeof(temp), "    -  %s : %.1f %% (%.1f%s of %.1f%s)", 
						counter->getName().c_str(), baseValue == 0 ? (Float) 0 : value/baseValue * 100, 
						value2, suffixesNumber[suffixIndex].c_str(),
						value3, suffixesNumber[suffixIndex].c_str());
					break;
				}
		}
		oss << temp << std::endl;
		++statsEntries;
	}

	if (statsEntries == 0) {
		oss << " * Statistics:" << std::endl
			<< "     none." << std::endl;
	}

	oss << "------------------------------------------------------------";
	m_mutex->unlock();
	return oss.str();
}

MTS_IMPLEMENT_CLASS(Statistics, false, Object)
MTS_NAMESPACE_END
