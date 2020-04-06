#version 330


in vec3 col;
in vec3 Position_worldspace;
in vec3 Normal_cameraspace;
in vec3 EyeDirection_cameraspace;
in vec3 LightPosition_worldspace;
in vec3 LightDirection_cameraspace;


out vec3 color;


uniform vec3 lightColor;
uniform float lightPower;

vec3 rgb2hsv(vec3 c)
{
	vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
	vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));
	vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));

	float d = q.x - min(q.w, q.y);
	float e = 1.0e-10;
	return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

vec3 hsv2rgb(vec3 c)
{
	vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
	vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
	return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

void main()
{


	vec3 rgb = hsv2rgb(col);


	vec3 MaterialDiffuseColor = rgb;
	vec3 MaterialAmbientColor = vec3(0.2, 0.2, 0.2) * MaterialDiffuseColor;
	vec3 MaterialSpecularColor = vec3(0.3, 0.3, 0.3);


	vec3 n = normalize(Normal_cameraspace);

	vec3 l = normalize(LightDirection_cameraspace);

	 float cosTheta = clamp(dot(n, l), 0, 1);

	vec3 E = normalize(EyeDirection_cameraspace);


	vec3 R = reflect(-l, n);

	 float cosAlpha = clamp(dot(E, R), 0, 1);

	color =

		MaterialAmbientColor +

		MaterialDiffuseColor * lightColor * cosTheta +

		MaterialSpecularColor * lightColor * pow(cosAlpha, 5);
}