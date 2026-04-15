// tinygltf implementation TU.
// stb_image is compiled in so tinygltf can decode embedded textures into
// the tinygltf::Image buffers that loadGLTF() surfaces on GltfScene::images.
// STB_IMAGE_STATIC keeps all symbols TU-local so we don't collide with any
// other stb_image TU the consuming app may already link (bro has its own).
#define STB_IMAGE_STATIC
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define TINYGLTF_IMPLEMENTATION
#define TINYGLTF_NO_INCLUDE_STB_IMAGE
#define TINYGLTF_NO_INCLUDE_STB_IMAGE_WRITE
#define TINYGLTF_NO_EXTERNAL_IMAGE
#define TINYGLTF_NO_STB_IMAGE_WRITE
#include "tiny_gltf.h"
