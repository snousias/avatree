#version 330



layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in vec3 vertexNormal_modelspace;
layout(location = 2) in vec3 vertexColor;


out vec3 col;
out vec3 Position_worldspace;
out vec3 Normal_cameraspace;
out vec3 EyeDirection_cameraspace;
out vec3 LightDirection_cameraspace;
out vec3 LightPosition_worldspace;


uniform mat4 projMatrix;
uniform mat4 mvMatrix;
uniform mat4 modelMatrix;
uniform mat4 viewMatrix;
uniform mat3 normalMatrix;
uniform float isWireframe;
uniform vec3 lightPosition_cameraspace;

void main()
{
	
	Position_worldspace = (modelMatrix * vec4(vertexPosition_modelspace, 1)).xyz;


	vec3 vertexPosition_cameraspace = (mvMatrix * vec4(vertexPosition_modelspace, 1)).xyz;
	EyeDirection_cameraspace = vec3(0, 0, 0) - vertexPosition_cameraspace;


	LightDirection_cameraspace = lightPosition_cameraspace + EyeDirection_cameraspace;
	LightPosition_worldspace = (inverse(viewMatrix) * vec4(lightPosition_cameraspace, 1)).xyz;

	
	Normal_cameraspace = normalMatrix * vertexNormal_modelspace;

	if (isWireframe > 0.5)
	{
		col = vec3(0.0, 0.0, 0.0);
		vec4 v = mvMatrix * vec4(vertexPosition_modelspace, 1);
		v.xyz *= 0.99;
	
		gl_Position = projMatrix * v;
	}
	else
	{
		col = vertexColor;

		gl_Position = projMatrix * mvMatrix * vec4(vertexPosition_modelspace, 1);
	}
}