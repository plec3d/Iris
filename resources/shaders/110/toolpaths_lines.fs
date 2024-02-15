#version 110

// Two key lights
const vec3 LIGHT_KEY_TOP_DIR = vec3(-1.0, 0.75, 1.0);
const vec3 LIGHT_KEY_BOT_DIR = vec3(-1.0, -0.75, 1.0);
#define LIGHT_KEY_DIFFUSE   0.16
#define LIGHT_KEY_SPECULAR  0.008

// One back light
const vec3 LIGHT_BACK_DIR = vec3(0.75, 0.5, -1.0);
#define LIGHT_BACK_DIFFUSE  0.3
#define LIGHT_BACK_SPECULAR  0.015

// Four fill lights
const vec3 LIGHT_FILL_TOP_FRONT_DIR = vec3(1.0, 0.75, 1.0);
const vec3 LIGHT_FILL_BOT_FRONT_DIR = vec3(1.0, -0.75, 1.0);
const vec3 LIGHT_FILL_TOP_MIDDLE_DIR = vec3(0.0, 1.0, 0.8);
const vec3 LIGHT_FILL_BOT_MIDDLE_DIR = vec3(0.0, -1.0, 0.8);
const vec3 LIGHT_FILL_CENTER_FRONT_DIR = vec3(0.0, 0.0, 1.0);

#define LIGHT_FILL_DIFFUSE  0.1
#define LIGHT_FILL_SPECULAR  0.004

#define LIGHT_SHININESS 7.0
#define INTENSITY_AMBIENT   0.3

// x = ambient, y = top diffuse, z = front diffuse, w = global
uniform vec4 light_intensity;
uniform vec4 uniform_color;

varying vec3 eye_normal;

void main()
{
    vec3 normal = normalize(eye_normal);

    // Compute the cos of the angle between the normal and lights direction. The light is directional so the direction is constant for every vertex.
    // Since these two are normalized the cosine is the dot product. Take the abs value to light the lines no matter in which direction the normal points.
    
    // Ambient light
    float intensity = INTENSITY_AMBIENT;

    // Key lights
    float NdotL = abs(dot(normal, LIGHT_KEY_TOP_DIR));
    intensity += NdotL * LIGHT_KEY_DIFFUSE;

    NdotL = abs(dot(normal, LIGHT_KEY_BOT_DIR));
    intensity += NdotL * LIGHT_KEY_DIFFUSE;

    // Back light
    NdotL = abs(dot(normal, LIGHT_BACK_DIR));
    intensity += NdotL * LIGHT_BACK_DIFFUSE;

    // Fill lights
    NdotL = abs(dot(normal, LIGHT_FILL_TOP_FRONT_DIR));
    intensity += NdotL * LIGHT_FILL_DIFFUSE;

    NdotL = abs(dot(normal, LIGHT_FILL_BOT_FRONT_DIR));
    intensity += NdotL * LIGHT_FILL_DIFFUSE;

    NdotL = abs(dot(normal, LIGHT_FILL_TOP_MIDDLE_DIR));
    intensity += NdotL * LIGHT_FILL_DIFFUSE;

    NdotL = abs(dot(normal, LIGHT_FILL_BOT_MIDDLE_DIR));
    intensity += NdotL * LIGHT_FILL_DIFFUSE;

    NdotL = abs(dot(normal, LIGHT_FILL_CENTER_FRONT_DIR));
    intensity += NdotL * LIGHT_FILL_DIFFUSE;

    gl_FragColor = vec4(uniform_color.rgb * 0.75 * intensity, uniform_color.a);
}
