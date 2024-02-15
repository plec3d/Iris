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

const vec3 ZERO = vec3(0.0, 0.0, 0.0);

uniform vec4 uniform_color;

varying vec3 clipping_planes_dots;
varying vec4 model_pos;

uniform bool volume_mirrored;

void main()
{
    if (any(lessThan(clipping_planes_dots, ZERO)))
        discard;
    vec3  color = uniform_color.rgb;
    float alpha = uniform_color.a;

    vec3 triangle_normal = normalize(cross(dFdx(model_pos.xyz), dFdy(model_pos.xyz)));
#ifdef FLIP_TRIANGLE_NORMALS
    triangle_normal = -triangle_normal;
#endif

    if (volume_mirrored)
        triangle_normal = -triangle_normal;

    // First transform the normal into camera space and normalize the result.
    vec3 eye_normal = normalize(gl_NormalMatrix * triangle_normal);
    // x = diffuse, y = specular
    vec2 intensity = vec2(0.0, 0.0);
    vec3 position = (gl_ModelViewMatrix * model_pos).xyz;

    // Compute the cos of the angle between the normal and lights direction. The light is directional so the direction is constant for every vertex.
    // Since these two are normalized the cosine is the dot product. We also need to clamp the result to the [0,1] range.
    
    // Ambient light
    intensity.x = INTENSITY_AMBIENT;

    // Key lights
    float NdotL = max(dot(eye_normal, LIGHT_KEY_TOP_DIR), 0.0);
    intensity.x += NdotL * LIGHT_KEY_DIFFUSE;
    intensity.y = LIGHT_KEY_SPECULAR * pow(max(dot(-normalize(position), reflect(-LIGHT_KEY_TOP_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    NdotL = max(dot(eye_normal, LIGHT_KEY_BOT_DIR), 0.0);
    intensity.x += NdotL * LIGHT_KEY_DIFFUSE;
    intensity.y += LIGHT_KEY_SPECULAR * pow(max(dot(-normalize(position), reflect(-LIGHT_KEY_BOT_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    // Back light
    NdotL = max(dot(eye_normal, LIGHT_BACK_DIR), 0.0);
    intensity.x += NdotL * LIGHT_BACK_DIFFUSE;
    intensity.y += LIGHT_BACK_SPECULAR * pow(max(dot(-normalize(position), reflect(-LIGHT_BACK_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    // Fill lights
    NdotL = max(dot(eye_normal, LIGHT_FILL_TOP_FRONT_DIR), 0.0);
    intensity.x += NdotL * LIGHT_FILL_DIFFUSE;
    intensity.y += LIGHT_FILL_SPECULAR * pow(max(dot(-normalize(position), reflect(-LIGHT_FILL_TOP_FRONT_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    NdotL = max(dot(eye_normal, LIGHT_FILL_BOT_FRONT_DIR), 0.0);
    intensity.x += NdotL * LIGHT_FILL_DIFFUSE;
    intensity.y += LIGHT_FILL_SPECULAR * pow(max(dot(-normalize(position), reflect(-LIGHT_FILL_BOT_FRONT_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    NdotL = max(dot(eye_normal, LIGHT_FILL_TOP_MIDDLE_DIR), 0.0);
    intensity.x += NdotL * LIGHT_FILL_DIFFUSE;
    intensity.y += LIGHT_FILL_SPECULAR * pow(max(dot(-normalize(position), reflect(-LIGHT_FILL_TOP_MIDDLE_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    NdotL = max(dot(eye_normal, LIGHT_FILL_TOP_MIDDLE_DIR), 0.0);
    intensity.x += NdotL * LIGHT_FILL_DIFFUSE;
    intensity.y += LIGHT_FILL_SPECULAR * pow(max(dot(-normalize(position), reflect(-LIGHT_FILL_TOP_MIDDLE_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    NdotL = max(dot(eye_normal, LIGHT_FILL_CENTER_FRONT_DIR), 0.0);
    intensity.x += NdotL * LIGHT_FILL_DIFFUSE;
    intensity.y += LIGHT_FILL_SPECULAR * pow(max(dot(-normalize(position), reflect(-LIGHT_FILL_CENTER_FRONT_DIR, eye_normal)), 0.0), LIGHT_SHININESS);

    gl_FragColor = vec4(vec3(intensity.y) + color * intensity.x, alpha);
}
