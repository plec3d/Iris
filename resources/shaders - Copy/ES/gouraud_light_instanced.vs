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

// vertex attributes
attribute vec3 v_position;
attribute vec3 v_normal;
// instance attributes
attribute vec3 i_offset;
attribute vec2 i_scales;

// x = diffuse, y = specular
varying vec2 intensity;

void main()
{
    // First transform the normal into camera space and normalize the result.
    vec3 eye_normal = normalize(gl_NormalMatrix * v_normal);
    vec4 world_position = vec4(v_position * vec3(vec2(1.5 * i_scales.x), 1.5 * i_scales.y) + i_offset - vec3(0.0, 0.0, 0.5 * i_scales.y), 1.0);
    vec3 eye_position = (gl_ModelViewMatrix * world_position).xyz;
    
    // Compute the cos of the angle between the normal and lights direction. The light is directional so the direction is constant for every vertex.
    // Since these two are normalized the cosine is the dot product. We also need to clamp the result to the [0,1] range.
    
    // Ambient light
    intensity.x = INTENSITY_AMBIENT;

    // Key lights
    float NdotL = max(dot(eye_normal, LIGHT_KEY_TOP_DIR), 0.0);
    intensity.x += NdotL * LIGHT_KEY_DIFFUSE;
    intensity.y = LIGHT_KEY_SPECULAR * pow(max(dot(-normalize(eye_position), reflect(-LIGHT_KEY_TOP_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    NdotL = max(dot(eye_normal, LIGHT_KEY_BOT_DIR), 0.0);
    intensity.x += NdotL * LIGHT_KEY_DIFFUSE;
    intensity.y += LIGHT_KEY_SPECULAR * pow(max(dot(-normalize(eye_position), reflect(-LIGHT_KEY_BOT_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    // Back light
    NdotL = max(dot(eye_normal, LIGHT_BACK_DIR), 0.0);
    intensity.x += NdotL * LIGHT_BACK_DIFFUSE;
    intensity.y += LIGHT_BACK_SPECULAR * pow(max(dot(-normalize(eye_position), reflect(-LIGHT_BACK_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    // Fill lights
    NdotL = max(dot(eye_normal, LIGHT_FILL_TOP_FRONT_DIR), 0.0);
    intensity.x += NdotL * LIGHT_FILL_DIFFUSE;
    intensity.y += LIGHT_FILL_SPECULAR * pow(max(dot(-normalize(eye_position), reflect(-LIGHT_FILL_TOP_FRONT_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    NdotL = max(dot(eye_normal, LIGHT_FILL_BOT_FRONT_DIR), 0.0);
    intensity.x += NdotL * LIGHT_FILL_DIFFUSE;
    intensity.y += LIGHT_FILL_SPECULAR * pow(max(dot(-normalize(eye_position), reflect(-LIGHT_FILL_BOT_FRONT_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    NdotL = max(dot(eye_normal, LIGHT_FILL_TOP_MIDDLE_DIR), 0.0);
    intensity.x += NdotL * LIGHT_FILL_DIFFUSE;
    intensity.y += LIGHT_FILL_SPECULAR * pow(max(dot(-normalize(eye_position), reflect(-LIGHT_FILL_TOP_MIDDLE_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    NdotL = max(dot(eye_normal, LIGHT_FILL_BOT_MIDDLE_DIR), 0.0);
    intensity.x += NdotL * LIGHT_FILL_DIFFUSE;
    intensity.y += LIGHT_FILL_SPECULAR * pow(max(dot(-normalize(eye_position), reflect(-LIGHT_FILL_BOT_MIDDLE_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    NdotL = max(dot(eye_normal, LIGHT_FILL_CENTER_FRONT_DIR), 0.0);
    intensity.x += NdotL * LIGHT_FILL_DIFFUSE;
    intensity.y += LIGHT_FILL_SPECULAR * pow(max(dot(-normalize(eye_position), reflect(-LIGHT_FILL_CENTER_FRONT_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    gl_Position = gl_ProjectionMatrix * vec4(eye_position, 1.0);
}
