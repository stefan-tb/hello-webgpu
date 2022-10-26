#define ZBUFFER
// see dawn/src/utils/ComboRenderPipelineDescriptor.cpp

#include "webgpu.h"
#include "matrix.h"

#include <string.h>
#include "mesh.h"

//float const data_positionAndNormals[] = {
//		-1.0f, -0.0f, 0.0f, 0.0f, 0.0f, 1.0f, // BL
//		 1.0f, -0.0f, 0.0f, 0.0f, 1.0f, 0.0f, // BR
//		-0.0f,  1.0f, 0.0f, 1.0f, 0.0f, 0.0f, // top
//};
//uint16_t const data_indices[] = {
//	0, 1, 2,
//	0 // padding (better way of doing this?)
//};

WGPUDevice device;
WGPUQueue queue;
WGPUSwapChain swapchain;

WGPURenderPipeline pipeline;
WGPUBuffer vertBuf; // vertex buffer with triangle position and colours
WGPUBuffer indxBuf; // index buffer
WGPUBuffer uniformBuf; // uniform buffer (containing the rotation angle)
WGPUBindGroup bindGroup;
#ifdef ZBUFFER
WGPUTexture depthStencilTexture;
#endif

/*
 * Workaround for Dawn currently expecting 'undefined' for entire buffers and
 * Emscripten/Chrome still expecting zero.
 */
#ifndef __EMSCRIPTEN__
#define ZERO_BUFFER_SIZE WGPU_WHOLE_SIZE
#else
#define ZERO_BUFFER_SIZE 0
#endif

 /**
  * Current rotation angle (in degrees, updated per frame).
  */
float rotDeg = 0.0f;

typedef struct Uniforms {
	matrix::Mat4x4 world;
	matrix::Mat4x4 viewProj;
} Uniforms;

/**
 * Vertex shader SPIR-V.
 * \code
 *	// glslc -Os -mfmt=num -o - -c in.vert
 *	#version 450
 *	layout(set = 0, binding = 0) uniform Rotation {
 *		float uRot;
 *	};
 *	layout(location = 0) in  vec2 aPos;
 *	layout(location = 1) in  vec3 aCol;
 *	layout(location = 0) out vec3 vCol;
 *	void main() {
 *		float cosA = cos(radians(uRot));
 *		float sinA = sin(radians(uRot));
 *		mat3 rot = mat3(cosA, sinA, 0.0,
 *					   -sinA, cosA, 0.0,
 *						0.0,  0.0,  1.0);
 *		gl_Position = vec4(rot * vec3(aPos, 1.0), 1.0);
 *		vCol = aCol;
 *	}
 * \endcode
 */
static uint32_t const triangle_vert_spirv[] = {
	0x07230203, 0x00010000, 0x000d0008, 0x00000043, 0x00000000, 0x00020011, 0x00000001, 0x0006000b,
	0x00000001, 0x4c534c47, 0x6474732e, 0x3035342e, 0x00000000, 0x0003000e, 0x00000000, 0x00000001,
	0x0009000f, 0x00000000, 0x00000004, 0x6e69616d, 0x00000000, 0x0000002d, 0x00000031, 0x0000003e,
	0x00000040, 0x00050048, 0x00000009, 0x00000000, 0x00000023, 0x00000000, 0x00030047, 0x00000009,
	0x00000002, 0x00040047, 0x0000000b, 0x00000022, 0x00000000, 0x00040047, 0x0000000b, 0x00000021,
	0x00000000, 0x00050048, 0x0000002b, 0x00000000, 0x0000000b, 0x00000000, 0x00050048, 0x0000002b,
	0x00000001, 0x0000000b, 0x00000001, 0x00050048, 0x0000002b, 0x00000002, 0x0000000b, 0x00000003,
	0x00050048, 0x0000002b, 0x00000003, 0x0000000b, 0x00000004, 0x00030047, 0x0000002b, 0x00000002,
	0x00040047, 0x00000031, 0x0000001e, 0x00000000, 0x00040047, 0x0000003e, 0x0000001e, 0x00000000,
	0x00040047, 0x00000040, 0x0000001e, 0x00000001, 0x00020013, 0x00000002, 0x00030021, 0x00000003,
	0x00000002, 0x00030016, 0x00000006, 0x00000020, 0x0003001e, 0x00000009, 0x00000006, 0x00040020,
	0x0000000a, 0x00000002, 0x00000009, 0x0004003b, 0x0000000a, 0x0000000b, 0x00000002, 0x00040015,
	0x0000000c, 0x00000020, 0x00000001, 0x0004002b, 0x0000000c, 0x0000000d, 0x00000000, 0x00040020,
	0x0000000e, 0x00000002, 0x00000006, 0x00040017, 0x00000018, 0x00000006, 0x00000003, 0x00040018,
	0x00000019, 0x00000018, 0x00000003, 0x0004002b, 0x00000006, 0x0000001e, 0x00000000, 0x0004002b,
	0x00000006, 0x00000022, 0x3f800000, 0x00040017, 0x00000027, 0x00000006, 0x00000004, 0x00040015,
	0x00000028, 0x00000020, 0x00000000, 0x0004002b, 0x00000028, 0x00000029, 0x00000001, 0x0004001c,
	0x0000002a, 0x00000006, 0x00000029, 0x0006001e, 0x0000002b, 0x00000027, 0x00000006, 0x0000002a,
	0x0000002a, 0x00040020, 0x0000002c, 0x00000003, 0x0000002b, 0x0004003b, 0x0000002c, 0x0000002d,
	0x00000003, 0x00040017, 0x0000002f, 0x00000006, 0x00000002, 0x00040020, 0x00000030, 0x00000001,
	0x0000002f, 0x0004003b, 0x00000030, 0x00000031, 0x00000001, 0x00040020, 0x0000003b, 0x00000003,
	0x00000027, 0x00040020, 0x0000003d, 0x00000003, 0x00000018, 0x0004003b, 0x0000003d, 0x0000003e,
	0x00000003, 0x00040020, 0x0000003f, 0x00000001, 0x00000018, 0x0004003b, 0x0000003f, 0x00000040,
	0x00000001, 0x0006002c, 0x00000018, 0x00000042, 0x0000001e, 0x0000001e, 0x00000022, 0x00050036,
	0x00000002, 0x00000004, 0x00000000, 0x00000003, 0x000200f8, 0x00000005, 0x00050041, 0x0000000e,
	0x0000000f, 0x0000000b, 0x0000000d, 0x0004003d, 0x00000006, 0x00000010, 0x0000000f, 0x0006000c,
	0x00000006, 0x00000011, 0x00000001, 0x0000000b, 0x00000010, 0x0006000c, 0x00000006, 0x00000012,
	0x00000001, 0x0000000e, 0x00000011, 0x0006000c, 0x00000006, 0x00000017, 0x00000001, 0x0000000d,
	0x00000011, 0x0004007f, 0x00000006, 0x00000020, 0x00000017, 0x00060050, 0x00000018, 0x00000023,
	0x00000012, 0x00000017, 0x0000001e, 0x00060050, 0x00000018, 0x00000024, 0x00000020, 0x00000012,
	0x0000001e, 0x00060050, 0x00000019, 0x00000026, 0x00000023, 0x00000024, 0x00000042, 0x0004003d,
	0x0000002f, 0x00000032, 0x00000031, 0x00050051, 0x00000006, 0x00000033, 0x00000032, 0x00000000,
	0x00050051, 0x00000006, 0x00000034, 0x00000032, 0x00000001, 0x00060050, 0x00000018, 0x00000035,
	0x00000033, 0x00000034, 0x00000022, 0x00050091, 0x00000018, 0x00000036, 0x00000026, 0x00000035,
	0x00050051, 0x00000006, 0x00000037, 0x00000036, 0x00000000, 0x00050051, 0x00000006, 0x00000038,
	0x00000036, 0x00000001, 0x00050051, 0x00000006, 0x00000039, 0x00000036, 0x00000002, 0x00070050,
	0x00000027, 0x0000003a, 0x00000037, 0x00000038, 0x00000039, 0x00000022, 0x00050041, 0x0000003b,
	0x0000003c, 0x0000002d, 0x0000000d, 0x0003003e, 0x0000003c, 0x0000003a, 0x0004003d, 0x00000018,
	0x00000041, 0x00000040, 0x0003003e, 0x0000003e, 0x00000041, 0x000100fd, 0x00010038
};

/**
 * WGSL equivalent of \c triangle_vert_spirv.
 */
static char const triangle_vert_wgsl[] = R"(
	let PI : f32 = 3.141592653589793;
	fn radians(degs : f32) -> f32 {
		return (degs * PI) / 180.0;
	}
	[[block]]
	struct VertexIn {
		[[location(0)]] Pos : vec3<f32>;
		[[location(1)]] Normal : vec3<f32>;
	};
	struct VertexOut {
		[[builtin(position)]] Position : vec4<f32>;
		[[location(0)]] Color : vec3<f32>;
	};
	[[block]]
	struct Uniforms {
		world : mat4x4<f32>;
		viewProj : mat4x4<f32>;
	};
	[[group(0), binding(0)]] var<uniform> uniforms : Uniforms;
	[[stage(vertex)]]
	fn main(input : VertexIn) -> VertexOut {
		let rot : mat4x4<f32> = mat4x4<f32>(
			vec4<f32>(1.0, 0.0, 0.0, 0.0),
			vec4<f32>(0.0, 1.0, 0.0, 0.0),
			vec4<f32>(0.0, 0.0, 1.0, 0.0),
			vec4<f32>(0.0, 0.0, 0.0, 1.0));

		var output : VertexOut;
		//output.Position = vec4<f32>(vec4<f32>(input.Pos, 1.0));
		output.Position = vec4<f32>(uniforms.world * vec4<f32>(input.Pos, 1.0));
		output.Position = vec4<f32>(uniforms.viewProj * output.Position);
		output.Color =  (uniforms.world * vec4<f32>(input.Normal, 0.0)).xyz;
		return output;
	}
)";

/**
 * Fragment shader SPIR-V.
 * \code
 *	// glslc -Os -mfmt=num -o - -c in.frag
 *	#version 450
 *	layout(location = 0) in  vec3 vCol;
 *	layout(location = 0) out vec4 fragColor;
 *	void main() {
 *		fragColor = vec4(vCol, 1.0);
 *	}
 * \endcode
 */
static uint32_t const triangle_frag_spirv[] = {
	0x07230203, 0x00010000, 0x000d0007, 0x00000013, 0x00000000, 0x00020011, 0x00000001, 0x0006000b,
	0x00000001, 0x4c534c47, 0x6474732e, 0x3035342e, 0x00000000, 0x0003000e, 0x00000000, 0x00000001,
	0x0007000f, 0x00000004, 0x00000004, 0x6e69616d, 0x00000000, 0x00000009, 0x0000000c, 0x00030010,
	0x00000004, 0x00000007, 0x00040047, 0x00000009, 0x0000001e, 0x00000000, 0x00040047, 0x0000000c,
	0x0000001e, 0x00000000, 0x00020013, 0x00000002, 0x00030021, 0x00000003, 0x00000002, 0x00030016,
	0x00000006, 0x00000020, 0x00040017, 0x00000007, 0x00000006, 0x00000004, 0x00040020, 0x00000008,
	0x00000003, 0x00000007, 0x0004003b, 0x00000008, 0x00000009, 0x00000003, 0x00040017, 0x0000000a,
	0x00000006, 0x00000003, 0x00040020, 0x0000000b, 0x00000001, 0x0000000a, 0x0004003b, 0x0000000b,
	0x0000000c, 0x00000001, 0x0004002b, 0x00000006, 0x0000000e, 0x3f800000, 0x00050036, 0x00000002,
	0x00000004, 0x00000000, 0x00000003, 0x000200f8, 0x00000005, 0x0004003d, 0x0000000a, 0x0000000d,
	0x0000000c, 0x00050051, 0x00000006, 0x0000000f, 0x0000000d, 0x00000000, 0x00050051, 0x00000006,
	0x00000010, 0x0000000d, 0x00000001, 0x00050051, 0x00000006, 0x00000011, 0x0000000d, 0x00000002,
	0x00070050, 0x00000007, 0x00000012, 0x0000000f, 0x00000010, 0x00000011, 0x0000000e, 0x0003003e,
	0x00000009, 0x00000012, 0x000100fd, 0x00010038
};

/**
 * WGSL equivalent of \c triangle_frag_spirv.
 */
static char const triangle_frag_wgsl[] = R"(
	[[stage(fragment)]]
	fn main([[location(0)]] normal : vec3<f32>) -> [[location(0)]] vec4<f32> {

	let surfaceNormal : vec3<f32> = normalize(normal);
	let color : vec3<f32> = vec3<f32>(1.0, 1.0, 1.0);
	let zero : vec3<f32> = vec3<f32>(0.0, 0.0, 0.0);
	let one : vec3<f32> = vec3<f32>(1.0, 1.0, 1.0);
	return vec4<f32>(
		color.xyz * 0.5 + // ambient (should be * 0.3)
		color.xyz * 0.4 * clamp(dot(surfaceNormal, vec3<f32>(0.0, 1.0, 0.0)), 0.0, 1.0) + // hemisperic
		color.xyz * 0.6 * clamp(dot(surfaceNormal, vec3<f32>(0.8660254, -0.4330127, 0.25)), 0.0, 1.0) * vec3<f32>(1.0, 0.9607844, 0.9078432)  // directional
		, 1.0);
	}
)";

/**
 * Helper to create a shader from SPIR-V IR.
 *
 * \param[in] code shader source (output using the \c -V \c -x options in \c glslangValidator)
 * \param[in] size size of \a code in bytes
 * \param[in] label optional shader name
 */
/*static*/ WGPUShaderModule createShader(const uint32_t* code, uint32_t size, const char* label = nullptr) {
	WGPUShaderModuleSPIRVDescriptor spirv = {};
	spirv.chain.sType = WGPUSType_ShaderModuleSPIRVDescriptor;
	spirv.codeSize = size / sizeof(uint32_t);
	spirv.code = code;
	WGPUShaderModuleDescriptor desc = {};
	desc.nextInChain = reinterpret_cast<WGPUChainedStruct*>(&spirv);
	desc.label = label;
	return wgpuDeviceCreateShaderModule(device, &desc);
}

/**
 * Helper to create a shader from WGSL source.
 *
 * \param[in] code WGSL shader source
 * \param[in] label optional shader name
 */
static WGPUShaderModule createShader(const char* const code, const char* label = nullptr) {
	WGPUShaderModuleWGSLDescriptor wgsl = {};
	wgsl.chain.sType = WGPUSType_ShaderModuleWGSLDescriptor;
	wgsl.source = code;
	WGPUShaderModuleDescriptor desc = {};
	desc.nextInChain = reinterpret_cast<WGPUChainedStruct*>(&wgsl);
	desc.label = label;
	return wgpuDeviceCreateShaderModule(device, &desc);
}

/**
 * Helper to create a buffer.
 *
 * \param[in] data pointer to the start of the raw data
 * \param[in] size number of bytes in \a data
 * \param[in] usage type of buffer
 */
static WGPUBuffer createBuffer(const void* data, size_t size, WGPUBufferUsage usage) {
	WGPUBufferDescriptor desc = {};
	desc.usage = WGPUBufferUsage_CopyDst | usage;
	desc.size = size;
	WGPUBuffer buffer = wgpuDeviceCreateBuffer(device, &desc);
	wgpuQueueWriteBuffer(queue, buffer, 0, data, size);
	return buffer;
}

/**
 * Bare minimum pipeline to draw a triangle using the above shaders.
 */
static void createPipelineAndBuffers() {
	// compile shaders
	// NOTE: these are now the WGSL shaders (tested with Dawn and Chrome Canary)
	WGPUShaderModule vertMod = createShader(triangle_vert_wgsl);
	WGPUShaderModule fragMod = createShader(triangle_frag_wgsl);

	// keep the old unused SPIR-V shaders around for a while...
	(void)triangle_vert_spirv;
	(void)triangle_frag_spirv;

	// bind group layout (used by both the pipeline layout and uniform bind group, released at the end of this function)
	WGPUBindGroupLayoutEntry bglEntry = {};
	bglEntry.binding = 0;
	bglEntry.visibility = WGPUShaderStage_Vertex;
	bglEntry.buffer.type = WGPUBufferBindingType_Uniform;

	WGPUBindGroupLayoutDescriptor bglDesc = {};
	bglDesc.entryCount = 1;
	bglDesc.entries = &bglEntry;
	WGPUBindGroupLayout bindGroupLayout = wgpuDeviceCreateBindGroupLayout(device, &bglDesc);

	// pipeline layout (used by the render pipeline, released after its creation)
	WGPUPipelineLayoutDescriptor layoutDesc = {};
	layoutDesc.bindGroupLayoutCount = 1;
	layoutDesc.bindGroupLayouts = &bindGroupLayout;
	WGPUPipelineLayout pipelineLayout = wgpuDeviceCreatePipelineLayout(device, &layoutDesc);

	// describe buffer layouts
	WGPUVertexAttribute vertAttrs[2] = {};
	vertAttrs[0].format = WGPUVertexFormat_Float32x3;
	vertAttrs[0].offset = 0;
	vertAttrs[0].shaderLocation = 0;
	vertAttrs[1].format = WGPUVertexFormat_Float32x3;
	vertAttrs[1].offset = 3 * sizeof(float);
	vertAttrs[1].shaderLocation = 1;
	WGPUVertexBufferLayout vertexBufferLayout = {};
	vertexBufferLayout.arrayStride = 6 * sizeof(float);
	vertexBufferLayout.attributeCount = 2;
	vertexBufferLayout.attributes = vertAttrs;

	// Fragment state
	WGPUBlendState blend = {};
	blend.color.operation = WGPUBlendOperation_Min;
	blend.color.srcFactor = WGPUBlendFactor_One;
	blend.color.dstFactor = WGPUBlendFactor_One;
	blend.alpha.operation = WGPUBlendOperation_Add;
	blend.alpha.srcFactor = WGPUBlendFactor_One;
	blend.alpha.dstFactor = WGPUBlendFactor_One;

	WGPUColorTargetState colorTarget = {};
	colorTarget.format = webgpu::getSwapChainFormat(device);
	colorTarget.blend = nullptr;//&blend;
	colorTarget.writeMask = WGPUColorWriteMask_All;

	// Set the defaults for the fragment state
	WGPUFragmentState fragment = {};
	fragment.module = fragMod;
	fragment.entryPoint = "main";
	fragment.targetCount = 1;
	fragment.targets = &colorTarget;

	WGPURenderPipelineDescriptor desc = {};
	desc.fragment = &fragment;

#ifdef ZBUFFER
	WGPUStencilFaceState stencilFace;
	stencilFace.compare = WGPUCompareFunction_Always;
	stencilFace.failOp = WGPUStencilOperation_Keep;
	stencilFace.depthFailOp = WGPUStencilOperation_Keep;
	stencilFace.passOp = WGPUStencilOperation_Keep;

	WGPUDepthStencilState depthStencilState = {};
	depthStencilState.format = WGPUTextureFormat_Depth24PlusStencil8;
	depthStencilState.depthWriteEnabled = true;
	depthStencilState.depthCompare = WGPUCompareFunction_Less;
	depthStencilState.stencilBack = stencilFace;
	depthStencilState.stencilFront = stencilFace;
	depthStencilState.stencilReadMask = 0xff;
	depthStencilState.stencilWriteMask = 0xff;
	depthStencilState.depthBias = 0;
	depthStencilState.depthBiasSlopeScale = 0.0;
	depthStencilState.depthBiasClamp = 0.0;
#endif

	// Other state
	desc.layout = pipelineLayout;
#ifdef ZBUFFER
	desc.depthStencil = &depthStencilState;
#endif

	// Set defaults for the vertex state.
	desc.vertex.module = vertMod;
	desc.vertex.entryPoint = "main";
	desc.vertex.bufferCount = 1;
	desc.vertex.buffers = &vertexBufferLayout;

	// Set the defaults for the multisample state
	desc.multisample.count = 1;
	desc.multisample.mask = 0xFFFFFFFF;
	desc.multisample.alphaToCoverageEnabled = false;

	// Set the defaults for the primitive state
	desc.primitive.frontFace = WGPUFrontFace_CCW;
	desc.primitive.cullMode = WGPUCullMode_None;
	desc.primitive.topology = WGPUPrimitiveTopology_TriangleList;
	desc.primitive.stripIndexFormat = WGPUIndexFormat_Undefined;

	pipeline = wgpuDeviceCreateRenderPipeline(device, &desc);

	// partial clean-up (just move to the end, no?)
	wgpuPipelineLayoutRelease(pipelineLayout);

	wgpuShaderModuleRelease(fragMod);
	wgpuShaderModuleRelease(vertMod);

	// create the buffers (x, y, r, g, b)
	vertBuf = createBuffer(data_positionAndNormals, sizeof(data_positionAndNormals), WGPUBufferUsage_Vertex);
	indxBuf = createBuffer(data_indices, sizeof(data_indices), WGPUBufferUsage_Index);

	// create the uniform bind group (note 'rotDeg' is copied here, not bound in any way)
	Uniforms dummy;
	uniformBuf = createBuffer(&dummy, sizeof(Uniforms), WGPUBufferUsage_Uniform);

	WGPUBindGroupEntry bgEntry = {};
	bgEntry.binding = 0;
	bgEntry.buffer = uniformBuf;
	bgEntry.offset = 0;
	bgEntry.size = sizeof(Uniforms);

	WGPUBindGroupDescriptor bgDesc = {};
	bgDesc.layout = bindGroupLayout;
	bgDesc.entryCount = 1;
	bgDesc.entries = &bgEntry;

	bindGroup = wgpuDeviceCreateBindGroup(device, &bgDesc);

	// last bit of clean-up
	wgpuBindGroupLayoutRelease(bindGroupLayout);

#ifdef ZBUFFER
	// https://github.com/samdauwe/webgpu-native-examples/blob/master/src/examples/msaa_line.c
	WGPUTextureDescriptor descriptor = {};
	descriptor.dimension = WGPUTextureDimension::WGPUTextureDimension_2D;
	descriptor.size.width = 800;
	descriptor.size.height = 450;
	descriptor.size.depthOrArrayLayers = 1;
	descriptor.sampleCount = 1;
	descriptor.format = WGPUTextureFormat::WGPUTextureFormat_Depth24PlusStencil8;
	descriptor.mipLevelCount = 1;
	descriptor.usage = WGPUTextureUsage::WGPUTextureUsage_RenderAttachment;
	depthStencilTexture = wgpuDeviceCreateTexture(device, &descriptor);
#endif
}

/**
 * Draws using the above pipeline and buffers.
 */
static bool redraw() {
	WGPUTextureView backBufView = wgpuSwapChainGetCurrentTextureView(swapchain);			// create textureView

#ifdef ZBUFFER
	//WGPUTextureViewDescriptor desc2{};
	//desc2.format = WGPUTextureFormat_Depth24PlusStencil8;
	//desc2.dimension = WGPUTextureViewDimension_2D;
	//desc2.baseMipLevel = 0;
	//desc2.mipLevelCount = 1;
	//desc2.baseArrayLayer = 0;
	//desc2.arrayLayerCount = 1;
	WGPUTextureView depthStencilView = wgpuTextureCreateView(depthStencilTexture, nullptr);
#endif

	WGPURenderPassColorAttachment colorAttachment = {};
	colorAttachment.loadOp = WGPULoadOp_Clear;
	colorAttachment.storeOp = WGPUStoreOp_Store;
	colorAttachment.clearColor.r = 0.3f;
	colorAttachment.clearColor.g = 0.3f;
	colorAttachment.clearColor.b = 0.3f;
	colorAttachment.clearColor.a = 1.0f;
	colorAttachment.view = backBufView;

#ifdef ZBUFFER
	WGPURenderPassDepthStencilAttachment depthStencilAttachment = {};
	depthStencilAttachment.clearDepth = 1.0f;
	depthStencilAttachment.clearStencil = 0;
	depthStencilAttachment.depthLoadOp = WGPULoadOp_Clear;
	depthStencilAttachment.depthStoreOp = WGPUStoreOp_Store;
	depthStencilAttachment.stencilLoadOp = WGPULoadOp_Clear;
	depthStencilAttachment.stencilStoreOp = WGPUStoreOp_Store;
	depthStencilAttachment.view = depthStencilView;
#endif

	WGPURenderPassDescriptor renderPass = {};
	renderPass.colorAttachmentCount = 1;
	renderPass.colorAttachments = &colorAttachment;
#ifdef ZBUFFER
	renderPass.depthStencilAttachment = &depthStencilAttachment;
#endif

	WGPUCommandEncoder encoder = wgpuDeviceCreateCommandEncoder(device, nullptr);			// create encoder
	WGPURenderPassEncoder pass = wgpuCommandEncoderBeginRenderPass(encoder, &renderPass);	// create pass

	// update the rotation, translate first, then rotate
	rotDeg += 0.5f;

	Uniforms uniforms;

	matrix::Mat4x4 m;
	matrix::Translate(uniforms.world, -25, 0, 0);
	matrix::Rotate(m, rotDeg, 1, 1, 0);
	matrix::Multiply(uniforms.world, m);
	matrix::Translate(m, 0, 0, -100);
	matrix::Multiply(uniforms.world, m);

	matrix::Ortho(uniforms.viewProj, 100 * (800 / 450.0), 100, -1000, 1000);

	matrix::Transpose(uniforms.world);
	matrix::Transpose(uniforms.viewProj);
	wgpuQueueWriteBuffer(queue, uniformBuf, 0, &uniforms, sizeof(uniforms));

	// draw the triangle (comment these five lines to simply clear the screen)
	wgpuRenderPassEncoderSetPipeline(pass, pipeline);
	wgpuRenderPassEncoderSetBindGroup(pass, 0, bindGroup, 0, 0);
	wgpuRenderPassEncoderSetVertexBuffer(pass, 0, vertBuf, 0, ZERO_BUFFER_SIZE);
	wgpuRenderPassEncoderSetIndexBuffer(pass, indxBuf, WGPUIndexFormat_Uint16, 0, ZERO_BUFFER_SIZE);
	wgpuRenderPassEncoderDrawIndexed(pass, (sizeof(data_indices) / sizeof(data_indices[0])), 1, 0, 0, 0);

	wgpuRenderPassEncoderEndPass(pass);
	wgpuRenderPassEncoderRelease(pass);														// release pass
	WGPUCommandBuffer commands = wgpuCommandEncoderFinish(encoder, nullptr);				// create commands
	wgpuCommandEncoderRelease(encoder);														// release encoder

	wgpuQueueSubmit(queue, 1, &commands);
	wgpuCommandBufferRelease(commands);														// release commands
#ifndef __EMSCRIPTEN__
	/*
	 * TODO: wgpuSwapChainPresent is unsupported in Emscripten, so what do we do?
	 */
	wgpuSwapChainPresent(swapchain);
#endif
#ifdef ZBUFFER
	wgpuTextureViewRelease(depthStencilView);												// release textureView
#endif
	wgpuTextureViewRelease(backBufView);													// release textureView

	return true;
}

extern "C" int __main__(int /*argc*/, char* /*argv*/[]) {
	if (window::Handle wHnd = window::create()) {
		if ((device = webgpu::create(wHnd))) {
			queue = wgpuDeviceGetQueue(device);
			swapchain = webgpu::createSwapChain(device);
			createPipelineAndBuffers();

			window::show(wHnd);
			window::loop(wHnd, redraw);

#ifndef __EMSCRIPTEN__
			wgpuBindGroupRelease(bindGroup);
			wgpuBufferRelease(uniformBuf);
			wgpuBufferRelease(indxBuf);
			wgpuBufferRelease(vertBuf);
			wgpuRenderPipelineRelease(pipeline);
			wgpuSwapChainRelease(swapchain);
			wgpuQueueRelease(queue);
			wgpuDeviceRelease(device);
#endif
		}
#ifndef __EMSCRIPTEN__
		window::destroy(wHnd);
#endif
	}
	return 0;
}
