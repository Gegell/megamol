/*
 * OpenGLException.cpp
 *
 * Copyright (C) 2006 by Universitaet Stuttgart (VIS). Alle Rechte vorbehalten.
 */

#ifdef _WIN32
#include <Windows.h>
#endif /* _WIN32 */

#include "vislib_gl/graphics/gl/OpenGLException.h"


#include <GL/glu.h>


/*
 * vislib_gl::graphics::gl::OpenGLException::OpenGLException
 */
vislib_gl::graphics::gl::OpenGLException::OpenGLException(const GLenum errorCode, const char* file, const int line)
        : Exception(file, line)
        , errorCode(errorCode) {
    Exception::setMsg(reinterpret_cast<const char*>(::gluErrorString(this->errorCode)));
}


/*
 * vislib_gl::graphics::gl::OpenGLException::OpenGLException
 */
vislib_gl::graphics::gl::OpenGLException::OpenGLException(const char* file, const int line)
        : Exception(file, line)
        , errorCode(::glGetError()) {
    Exception::setMsg(reinterpret_cast<const char*>(::gluErrorString(this->errorCode)));
}


/*
 * vislib_gl::graphics::gl::OpenGLException::OpenGLException
 */
vislib_gl::graphics::gl::OpenGLException::OpenGLException(const OpenGLException& rhs)
        : Exception(rhs)
        , errorCode(rhs.errorCode) {}


/*
 * vislib_gl::graphics::gl::OpenGLException::~OpenGLException
 */
vislib_gl::graphics::gl::OpenGLException::~OpenGLException(void) {}


/*
 * vislib_gl::graphics::gl::OpenGLException::operator =
 */
vislib_gl::graphics::gl::OpenGLException& vislib_gl::graphics::gl::OpenGLException::operator=(
    const OpenGLException& rhs) {
    Exception::operator=(rhs);

    if (this != &rhs) {
        this->errorCode = rhs.errorCode;
    }

    return *this;
}
