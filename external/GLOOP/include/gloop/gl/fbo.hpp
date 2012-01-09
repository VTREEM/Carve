#pragma once

#include <vector>
#include <string>

namespace gloop {
  class FrameBuffer {
    static std::vector<GLuint> s_bound_fbo;
    
    GLuint fbo;
    unsigned width;
    unsigned height;
    
    std::pair<GLenum, Surface::Ptr> depth;
    std::pair<GLenum, Surface::Ptr> stencil;
    std::vector<std::pair<GLenum, Surface::Ptr> > colour;
    
    bool _attach(GLenum target, const Surface::Ptr &surf, GLenum attachment, int mipmap_level);
    
  public:
      FrameBuffer(unsigned _width, unsigned _height);
    
    void add(GLenum target, const Surface::Ptr &surf);
    void add(const Surface::Ptr &surf);
    void init();
    void pushBind();
    void bind();
    static void popBind();
    void attach(int mipmap_level = 0);
  };
}