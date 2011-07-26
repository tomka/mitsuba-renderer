
/**********************************************************************\
* AUTHOR : HILLAIRE Sébastien
*
* MAIL   : hillaire_sebastien@yahoo.fr
* SITE   : sebastien.hillaire.free.fr
*
*	You are free to totally or partially use this file/code.
* If you do, please credit me in your software or demo and leave this
* note.
*	Share your work and your ideas as much as possible!
\*********************************************************************/

#include <mitsuba/mitsuba.h>
#if defined(__OSX__)
#include <OpenGL/glew.h>
#else
#include <GL/glew.h>
#endif
#include <mitsuba/hw/framebufferobject.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>

#ifndef GL_FRAMEBUFFER_INCOMPLETE_DUPLICATE_ATTACHMENT_EXT
#define GL_FRAMEBUFFER_INCOMPLETE_DUPLICATE_ATTACHMENT_EXT 0x8CD8
#endif

MTS_NAMESPACE_BEGIN

bool CheckFramebufferStatus(const GLuint fbo)
{
    GLenum status;
	bool ret = false;

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);

    status = (GLenum) glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    switch(status) {
        case GL_FRAMEBUFFER_COMPLETE_EXT:
			std::cerr << "Frame Buffer Created" << std::endl;
			ret = true;
            break;
        case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
			std::cerr << "Unsupported framebuffer format" << std::endl;
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
            std::cerr << "Framebuffer incomplete, missing attachment" << std::endl;
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_DUPLICATE_ATTACHMENT_EXT:
            std::cerr << "Framebuffer incomplete, duplicate attachment" << std::endl;
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
            std::cerr << "Framebuffer incomplete, attached images must have same dimensions" << std::endl;
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
            std::cerr << "Framebuffer incomplete, attached images must have same format" << std::endl;
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
            std::cerr << "Framebuffer incomplete, missing draw buffer" << std::endl;
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
            std::cerr << "Framebuffer incomplete, missing read buffer" << std::endl;
            break;
        default:
            std::cerr << "Framebuffer incomplete, UNKNOW ERROR" << std::endl;
            //Assert(false);
    }

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	return ret;
}

//////////////////////////////////////////////////////////////////////////////////

FrameBufferObject::FrameBufferObject(unsigned int nbColorBuffer)
        : nbColorAttachement(nbColorBuffer)
{
	if(nbColorAttachement > getMaxColorAttachments())
		nbColorAttachement = getMaxColorAttachments();
}

void FrameBufferObject::init(unsigned int width, unsigned int height,
		const unsigned int* colorBufferInternalFormat, 
		const unsigned int* colorBufferSWRAP,
		const unsigned int* colorBufferTWRAP,
		const unsigned int* colorBufferMinFiltering,
		const unsigned int* colorBufferMagFiltering,
		FBO_DepthBufferType depthBufferType,
		const unsigned int depthBufferMinFiltering,
		const unsigned int depthBufferMagFiltering,
		const unsigned int depthBufferSWRAP,
		const unsigned int depthBufferTWRAP,
		bool depthTextureCompareToR) {
	
	/////////////////INITIALIZATION/////////////////
    this->width = width;
    this->height = height;

	//color render buffer
	if(this->nbColorAttachement>0)
	{
		colorTextures = new GLuint[nbColorAttachement];
		colorMinificationFiltering = new GLuint[nbColorAttachement];
		for(int i=0; i<this->nbColorAttachement; i++)
		{
			glGenTextures(1, &colorTextures[i]);
			glBindTexture(GL_TEXTURE_2D, colorTextures[i]);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, colorBufferMinFiltering[i]);
			this->colorMinificationFiltering[i] = colorBufferMinFiltering[i];
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, colorBufferMagFiltering[i]);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, colorBufferSWRAP[i]);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, colorBufferTWRAP[i]);
			glTexImage2D(GL_TEXTURE_2D, 0, colorBufferInternalFormat[i], width, height, 0, GL_RGBA, GL_FLOAT, NULL);

			if(this->colorMinificationFiltering[i]==GL_NEAREST_MIPMAP_NEAREST
			|| this->colorMinificationFiltering[i]==GL_LINEAR_MIPMAP_NEAREST
			||this->colorMinificationFiltering[i]==GL_NEAREST_MIPMAP_LINEAR
			|| this->colorMinificationFiltering[i]==GL_LINEAR_MIPMAP_LINEAR)
			{
				glGenerateMipmapEXT(GL_TEXTURE_2D);
			}
		}
	}

	//depth render buffer
	this->depthType = depthBufferType;
	if(this->depthType!=FBO_DepthBufferType_NONE)
	{
		switch(this->depthType)
		{
		case FBO_DepthBufferType_TEXTURE:
			glGenTextures(1, &this->depthID);
			glBindTexture(GL_TEXTURE_2D, this->depthID);
			this->depthMinificationFiltering = depthBufferMinFiltering;
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, depthBufferMinFiltering);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, depthBufferMagFiltering);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, depthBufferSWRAP);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, depthBufferTWRAP);
			if(depthTextureCompareToR)
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE_ARB, GL_COMPARE_R_TO_TEXTURE_ARB);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24_ARB, width, height, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, NULL);

			if(this->depthMinificationFiltering==GL_NEAREST_MIPMAP_NEAREST
			|| this->depthMinificationFiltering==GL_LINEAR_MIPMAP_NEAREST
			||this->depthMinificationFiltering==GL_NEAREST_MIPMAP_LINEAR
			|| this->depthMinificationFiltering==GL_LINEAR_MIPMAP_LINEAR)
			{
				glGenerateMipmapEXT(GL_TEXTURE_2D);
			}
			break;

		case FBO_DepthBufferType_RENDERTARGET:
		default:
			glGenRenderbuffersEXT(1, &this->depthID);
			glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, this->depthID);
			glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT24_ARB, width, height);
			break;
		}
	}

	/////////////////ATTACHEMENT/////////////////
    glGenFramebuffersEXT(1, &fbo);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);

	//color render buffer attachement
	if(nbColorAttachement>0)
	{
		for(int i=0; i<this->nbColorAttachement; i++)
		{
			glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT+i, GL_TEXTURE_2D, this->colorTextures[i], 0 );
		}
		glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
		glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
	}
	else
	{
		glReadBuffer(GL_NONE);
		glDrawBuffer(GL_NONE);
	}
	
	//depth render buffer attachement
	if(this->depthType!=FBO_DepthBufferType_NONE)
	{
		switch(this->depthType)
		{
		case FBO_DepthBufferType_TEXTURE:
			glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, this->depthID, 0);
			break;

		case FBO_DepthBufferType_RENDERTARGET:
		default:
			glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, this->depthID);
			break;
		}
	}

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	CheckFramebufferStatus(this->fbo);
}

FrameBufferObject::~FrameBufferObject()
{
	if(this->nbColorAttachement>0)
	{
		delete [] this->colorTextures;
		delete [] this->colorMinificationFiltering;
	}
}



//////////////////////////////////////////////////////////////////////////////////



const GLint FrameBufferObject::getMaxColorAttachments()
{
  GLint maxAttach = 0;
  glGetIntegerv( GL_MAX_COLOR_ATTACHMENTS_EXT, &maxAttach );
  return maxAttach;
}
const GLint FrameBufferObject::getMaxBufferSize()
{
  GLint maxSize = 0;
  glGetIntegerv( GL_MAX_RENDERBUFFER_SIZE_EXT, &maxSize );
  return maxSize;
}



void FrameBufferObject::enableRenderToColorAndDepth(const unsigned int colorBufferNum) const
{
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);

	if (nbColorAttachement > 0)
	{
		unsigned int colorBuffer = GL_COLOR_ATTACHMENT0_EXT + colorBufferNum;
		if( (int)colorBufferNum > nbColorAttachement)
			colorBuffer = GL_COLOR_ATTACHMENT0_EXT;
		glDrawBuffer(colorBuffer);
	}
	else
	{
		glDrawBuffer(GL_NONE);	//for example, in the case of rendering in a depth buffer
	}
}

void FrameBufferObject::enableRenderToColorAndDepth_MRT(const GLuint numBuffers, const GLenum* drawbuffers) const
{ 
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, this->fbo);

    glDrawBuffers(numBuffers, drawbuffers);
	//CheckFramebufferStatus(this->fbo);
}

void FrameBufferObject::disableRenderToColorDepth()
{
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}



void FrameBufferObject::saveAndSetViewPort() const
{
	glPushAttrib(GL_VIEWPORT_BIT);
    glViewport(0, 0, this->width, this->height);
}

void FrameBufferObject::restoreViewPort() const
{
	glPopAttrib();
}



void FrameBufferObject::bindColorTexture(const unsigned int colorBufferNum) const
{
	if (nbColorAttachement > 0 && (int)colorBufferNum < nbColorAttachement) {
		glBindTexture(GL_TEXTURE_2D, colorTextures[colorBufferNum]);
	} else {
		glBindTexture(GL_TEXTURE_2D, 0);
    }
}

void FrameBufferObject::bindDepthTexture() const
{
	if(this->depthType==FBO_DepthBufferType_TEXTURE)
	{
		glBindTexture(GL_TEXTURE_2D, this->depthID);
	}
	else
		glBindTexture(GL_TEXTURE_2D, 0);
}

void FrameBufferObject::generateColorBufferMipMap(const unsigned int colorBufferNum) const
{
	if(this->nbColorAttachement>0 && (int)colorBufferNum<this->nbColorAttachement)
	{
		if(this->colorMinificationFiltering[colorBufferNum]==GL_NEAREST
		|| this->colorMinificationFiltering[colorBufferNum]==GL_LINEAR)
			return;	//don't allow to generate mipmap chain for texture that don't support it at the creation

		glBindTexture(GL_TEXTURE_2D, this->colorTextures[colorBufferNum]);
		glGenerateMipmapEXT(GL_TEXTURE_2D);
	}
}

void FrameBufferObject::generateDepthBufferMipMap() const
{
	if(this->depthType==FBO_DepthBufferType_TEXTURE)
	{
		if(this->depthMinificationFiltering==GL_NEAREST
		|| this->depthMinificationFiltering==GL_LINEAR)
			return;	//don't allow to generate mipmap chain for texture that don't support it at the creation

		glBindTexture(GL_TEXTURE_2D, this->depthID);
		glGenerateMipmapEXT(GL_TEXTURE_2D);
	}
}



unsigned int FrameBufferObject::getNumberOfColorAttachement() const
{
	return this->nbColorAttachement;
}

FBO_DepthBufferType FrameBufferObject::getDepthBufferType() const
{
	return this->depthType;
}

unsigned int FrameBufferObject::getWidth() const
{
	return this->width;
}

unsigned int FrameBufferObject::getHeight() const
{
	return this->height;
}

void FrameBufferObject::saveToDisk(const unsigned int colorBufferNum, const std::string &path) const {
	if(this->nbColorAttachement>0 && (int)colorBufferNum<this->nbColorAttachement)
	{
        float *data = new float[width * height * 3];
        glEnable(GL_TEXTURE_2D);
        bindColorTexture(colorBufferNum);
        glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_FLOAT, data);
        glDisable(GL_TEXTURE_2D);

        ref<Bitmap> bmp = new Bitmap(width, height, 128);
            for(unsigned int i=0; i < width*height; i+=1) {
                bmp->getFloatData()[i*4] = data[i*3];
                bmp->getFloatData()[i*4+1] = data[i*3+1];
                bmp->getFloatData()[i*4+2] = data[i*3+2];
                bmp->getFloatData()[i*4+3] = 1.0f;
            }
         bmp->save(Bitmap::EEXR, new FileStream(path, FileStream::ETruncWrite) );
        delete [] data;
	}
}

MTS_IMPLEMENT_CLASS(FrameBufferObject, false, Object)
MTS_NAMESPACE_END
