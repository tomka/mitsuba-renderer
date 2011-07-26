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

#ifndef QTGUI_COMMON_H
#define QTGUI_COMMON_H

#include <mitsuba/core/platform.h>
#include <QtGui>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/vpl.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/bitmap.h>
#include "model/snowproperties.h"
#include "snowmaterialmanager.h"

using namespace mitsuba;

enum EConnectionType {
	ESSHConnection = 0,
	EDirectConnection
};

enum ENavigationMode {
	EStandard = 0,
	EFlythrough
};

enum ESelectionMode {
	ENothing = 0,
	EShape,
	EScene
};

enum EGeneralRenderMode {
    ERealtime = 0,
    EOffline
};

/* Different surface rendering modes */
enum ESurfaceRenderMode {
    ENoSurface = 0,
    EWiscombeWarrenAlbedo,
    EWiscombeWarrenBRDF,
    EHanrahanKruegerBRDF,
    EMicrofacetBRDF
};

/* Different sub-surface rendering modes */
enum ESubSurfaceRenderMode {
    ENoSubSurface = 0,
    EJensenDipoleBSSRDF,
    EJensenMultipoleBSSRDF,
    EJakobADipoleBSSRDF
};

Q_DECLARE_METATYPE( Shape * );

struct SnowRenderSettings {
    EGeneralRenderMode generalRenderMode;
    ESurfaceRenderMode surfaceRenderMode;
    ESubSurfaceRenderMode subsurfaceRenderMode;

    /* Wiscombe BRDF settings */
    Float wiscombeDepth;

    /* Hanrahan-Krueger BRDF settings */
    Float hkSingleScatteringFactor;
    Float hkMultipleScatteringFactor;
    bool hkUseMultipleScattering;

    /* Jensen dipole BSSRDF settings */
    Float dipoleDensityFactor;
    Float dipoleSampleFactor;
    bool dipoleUseSingleScattering;
    bool dipoleMartelliDC;
    bool dipoleTexture;
    std::string dipoleZrTexture;
    std::string dipoleSigmaTrTexture;
    Float dipoleTextureUScaling;
    Float dipoleTextureVScaling;
    bool dipoleDumpIrrtree;
    std::string dipoleDumpIrrtreePath;
    bool dipoleUseLut;
    Float dipoleLutResolution;
    int dipoleLutMCIterations;
    Float dipoleLutRmax;
    bool dipoleLutPredefineRmax;
    bool dipoleHasRoughSurface;

    /* Jensen multipole settings */
    Float multipoleDensityFactor;
    Float multipoleSampleFactor;
    int multipoleExtraDipoles;
    Float multipoleSlabThickness;
    bool multipoleUseSingleScattering;
    bool multipoleMartelliDC;
    bool multipoleUseLut;
    Float multipoleLutResolution;
    int multipoleLutMCIterations;
    Float multipoleLutRmax;
    bool multipoleLutPredefineRmax;

    /* Jakob anisotropic dipole settings */
    Float adipoleDensityFactor;
    Float adipoleSampleFactor;
    Float adipoleSigmaTn;
    std::string adipoleD;

    /* Shah realtime SSS */
    enum EShahAlbedoType { EWiscombeWarrenAlbedo = 0, EWhiteAlbedo = 1, ECustomAlbedo = 2};
    enum EShahDiffusionPrType { ESnowProfile = 0, EExampleProfile = 1};

    bool shahExpandSilhouette;
    bool shahShowSplatOrigins;
    bool shahShowLight;
    EShahAlbedoType shahAlbedoMapType;
    std::string shahAlbedoMapCustomPath;
    ref<Bitmap> shahAlbedoMap;
    EShahDiffusionPrType shahDiffusionProfileType;
    int shahDiffusionExample;
    ref<Bitmap> shahDiffusionProfile;
    Float shahRmax;
    int shahMCIterations;
    bool shahPredefineRmax;
    Spectrum shahSpecularColor;
    Float shahErrorThreshold;
    int shahMaxLightViewResolution; 
    int shahBackbufferWidth, shahBackbufferHeight;
    Float shahWeight;
    Float shahExposure;

    SnowRenderSettings() :
        generalRenderMode(EOffline), surfaceRenderMode(ENoSurface), subsurfaceRenderMode(ENoSubSurface),
        wiscombeDepth(2.0f), hkSingleScatteringFactor(1.0f), hkMultipleScatteringFactor(1.0f),
        dipoleDensityFactor(1.0f), dipoleSampleFactor(6.0f), dipoleUseSingleScattering(false),
        dipoleMartelliDC(false), dipoleTexture(false), dipoleTextureUScaling(1.0f), dipoleTextureVScaling(1.0f),
        dipoleDumpIrrtree(false), dipoleDumpIrrtreePath(""),
        dipoleUseLut(true), dipoleLutResolution(0.001), dipoleLutMCIterations(10000), dipoleLutRmax(10.0f),
        dipoleLutPredefineRmax(true), dipoleHasRoughSurface(false),
        multipoleDensityFactor(1.0f), multipoleSampleFactor(6.0f), multipoleExtraDipoles(2),
        multipoleSlabThickness(0.2f), multipoleUseSingleScattering(false), multipoleMartelliDC(false),
        multipoleUseLut(true), multipoleLutResolution(0.001), multipoleLutMCIterations(10000),
        multipoleLutRmax(10.0f), multipoleLutPredefineRmax(true),
        adipoleDensityFactor(1.0f), adipoleSampleFactor(1.0f), adipoleSigmaTn(1.0f),
        // default to sin^20 flake distribution
        adipoleD("0.47827, 7.5057e-09, -4.313e-10, 7.5057e-09, 0.47827, 2.5069e-10, -4.313e-10, 2.5069e-10, 0.043454"),
        shahExpandSilhouette(true), shahShowSplatOrigins(false), shahShowLight(false),
        shahAlbedoMapType(EWiscombeWarrenAlbedo), shahDiffusionProfileType(ESnowProfile), shahDiffusionExample(0),
        shahRmax(0.5f), shahMCIterations(10000), shahPredefineRmax(true), shahSpecularColor(0.0f),
        shahErrorThreshold(0.01), shahMaxLightViewResolution(1024), shahBackbufferWidth(100),
        shahBackbufferHeight(75), shahWeight(1.0f), shahExposure(0.1)
    {
        /* try to load last texture paths */
	    QSettings settings("mitsuba-renderer.org", "qtgui");
		dipoleZrTexture = settings.value("lastDipoleZrTexture", "").toString().toStdString();
		dipoleSigmaTrTexture = settings.value("lastDipoleSigmaTrTexture", "").toString().toStdString();
        /* load default albedo and diffusion maps for realtime SSS */
        {   QResource res("/resources/snow/white.bmp");
            SAssert(res.isValid());
            ref<Stream> mStream = new MemoryStream(res.size());
            mStream->write(res.data(), res.size());
            mStream->setPos(0);
            ref<Bitmap> bitmap = new Bitmap(Bitmap::EBMP, mStream);
            shahAlbedoMap = bitmap; }
        {   std::ostringstream name; name << "/resources/snow/diffProfExample4.bmp";
            QResource res(name.str().c_str());
            SAssert(res.isValid());
            ref<Stream> mStream = new MemoryStream(res.size());
            mStream->write(res.data(), res.size());
            mStream->setPos(0);
            ref<Bitmap> bitmap = new Bitmap(Bitmap::EBMP, mStream);
            shahDiffusionProfile = bitmap; }
    }
};

namespace mitsuba  {
	class RemoteWorker;
};

struct ServerConnection {
	EConnectionType type;
	QString hostName, userName, instDir;
	int port;
	RemoteWorker *worker;
	bool isRegistered;

	inline ServerConnection() : worker(NULL), isRegistered(false) { }

	inline bool operator==(const ServerConnection &c) const {
		return type == c.type && hostName == c.hostName 
			&& userName == c.userName && instDir == c.instDir
			&& port == c.port && worker == c.worker
			&& isRegistered == c.isRegistered;
	}

	inline void fromVariant(QList<QVariant> list) {
		type = (EConnectionType) list[0].toInt();
		hostName = list[1].toString();
		port = list[2].toInt();
		if (type == ESSHConnection) {
			userName = list[3].toString();
			instDir = list[4].toString();
		}
	}

	inline QList<QVariant> toVariant() const {
		QList<QVariant> result;
		result.append(type);
		result.append(hostName);
		result.append(port);
		if (type == ESSHConnection) {
			result.append(userName);
			result.append(instDir);
		}
		return result;
	}

	inline QByteArray toByteArray() const {
		QByteArray a;
		QDataStream stream(&a, QIODevice::WriteOnly);
		stream << toVariant();
		return a;
	}
	
	inline void fromByteArray(QByteArray a) {
		QDataStream stream(a);
		QList<QVariant> variant;
		stream >> variant;
		fromVariant(variant);
	}

	bool createWorker(QWidget *parent);
	
	QString toString() const;
};

enum EMode {
	EPreview = 0,
	ERender
};

enum EPreviewMethod {
	EDisabled = 0,
	EOpenGL,
	EOpenGLSinglePass,
	ERayTraceCoherent,
	ERayTrace,
    EOpenGLRealtime
};

enum EToneMappingMethod {
	EGamma = 0,
	EReinhard
};

namespace mitsuba {
	class GPUSync;
	class GPUTexture;
};

struct PreviewQueueEntry {
	int id;
	size_t vplSampleOffset;
	GPUTexture *buffer;
	GPUSync *sync;

	inline PreviewQueueEntry(int id = 0) 
		: id(id), vplSampleOffset(0), buffer(NULL), sync(NULL) {
	}
};

struct VisualWorkUnit {
	Point2i offset;
	Vector2i size;
	int worker;
};

struct block_comparator : std::binary_function<VisualWorkUnit, VisualWorkUnit, bool> {
	static int compare(const VisualWorkUnit &v1, const VisualWorkUnit &v2) {
		if (v1.offset.x < v2.offset.x) return -1;
		else if (v1.offset.x > v2.offset.x) return 1;
		if (v1.offset.y < v2.offset.y) return -1;
		else if (v1.offset.y > v2.offset.y) return 1;
		if (v1.size.x < v2.size.x) return -1;
		else if (v1.size.x > v2.size.x) return 1;
		if (v1.size.y < v2.size.y) return -1;
		else if (v1.size.y > v2.size.y) return 1;
		return 0;
	}

	bool operator()(const VisualWorkUnit &v1, const VisualWorkUnit &v2) const {
		return compare(v1, v2) < 0;
	}
};

struct SceneContext {
	/* Scene-related */
	ref<Scene> scene;
	int sceneResID;
	QString fileName;
	QString shortName;
	Float movementScale;
	Vector up;

	/* Rendering/Preview-related */
	RenderJob *renderJob;
	bool cancelled;
	float progress;
	QString eta, progressName;
	ref<Bitmap> framebuffer;
	std::set<VisualWorkUnit, block_comparator> workUnits;
	EMode mode, cancelMode;
	Float gamma, exposure, clamping;
	bool srgb;
	int pathLength, shadowMapResolution;
	EPreviewMethod previewMethod;
	EToneMappingMethod toneMappingMethod;
	QSize windowSize, sizeIncrease;
	Vector2i scrollOffset;
	Float reinhardKey, reinhardBurn;
	bool diffuseSources;
	bool diffuseReceivers;
	bool showKDTree;
	int shownKDTreeLevel;
	ESelectionMode selectionMode;
	const Shape *selectedShape;
	bool showNormals;
    Float normalScaling;

	/* Preview state */
	std::deque<VPL> vpls;
	PreviewQueueEntry previewBuffer;

	SceneContext() : scene(NULL), sceneResID(-1), 
		renderJob(NULL), selectionMode(ENothing),
		selectedShape(NULL),
        normalScaling(0.04), currentlySelectedShape(NULL)
    { }

    /* Snow properties */
    SnowProperties snow;
    SnowRenderSettings snowRenderSettings;
    SnowMaterialManager snowMaterialManager;

    /* the currently selected shape */
    Shape *currentlySelectedShape;

	/// Detect the path length
	int detectPathLength() const;

	/* Clone a scene */
	SceneContext(SceneContext *ctx);
	~SceneContext();
};

class NonClosableDialog : public QDialog {
public:
	NonClosableDialog(QWidget *parent) : QDialog(parent) {
	}

	void closeEvent(QCloseEvent *e) {
		e->ignore();
	}
};

class ProgramVersion {
public:
	inline ProgramVersion(const QString &versionString) {
		QStringList sl;
		sl = versionString.trimmed().split('.');
		SAssert(sl.size() == 3);
		major = sl[0].toInt();
		minor = sl[1].toInt();
		release = sl[2].toInt();
	}

	inline bool operator<(const ProgramVersion &other) const {
		if (major < other.major)
			return true;
		if (major > other.major)
			return false;
		if (minor < other.minor)
			return true;
		if (minor > other.minor)
			return false;
		if (release < other.release)
			return true;
		return false;
	}

	inline bool operator==(const ProgramVersion &other) const {
		return major == other.major 
			&& minor == other.minor 
			&& release == other.release;
	}

	QString toString() const {
		return QString("%1.%2.%3").arg(major).arg(minor).arg(release);
	}
private:
	int major;
	int minor;
	int release;
};

#endif // QTGUI_COMMON_H
