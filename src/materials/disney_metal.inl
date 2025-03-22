#include "../microfacet.h"
#include "matrix.h"


Real GTR(Spectrum v, Real ax, Real ay) {
  Real A =
      (sqrt(1.0 + ((pow(v.x * ax, 2.0) + pow(v.y * ay, 2.0)) / pow(v.z, 2.0))) -
       1.0) /
      2.0;
  return 1.0 / (1.0 + A);
}

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    Real n_dot_out = dot(frame.n, dir_out);
    if (n_dot_out <= 0 ) {
        return make_zero_spectrum();
    }
    auto baseColor =
        eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    auto roughness =
        eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    auto anisotropic =
        eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    auto h = normalize(dir_in + dir_out);
    Spectrum Fresnel =
        baseColor + (1.0 - baseColor) * pow((1.0 - dot(h, dir_out)), 5.0);
    Real aspect = sqrt(1.0 - anisotropic * 0.9);
    Real ax = max(0.001, roughness * roughness / aspect);
    Real ay = max(0.001, roughness * roughness * aspect);
    // Real G = GTR(dir_in, ax, ay) * GTR(dir_out, ax, ay);
    Real G = GTR(dir_in, ax, ay) * GTR(dir_out, ax, ay);
    Real hlx = dot(frame.x, h);
    Real hly = dot(frame.y, h);
    Real hlz = dot(frame.n, h);
    Real hlxyz = ((hlx * hlx) / pow(ax, 2.0) + (hlx * hlx) / pow(ax, 2.0) +
                  (hlx * hlx) / pow(ax, 2.0)) *
                 ((hlx * hlx) / pow(ax, 2.0) + (hlx * hlx) / pow(ax, 2.0) +
                  (hlx * hlx) / pow(ax, 2.0));
    Real D = 1.0 / (c_PI * ax * ay * hlxyz);
    return Fresnel * G * D / (4.0 * dot(vertex.shading_frame.n, dir_in));
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    auto baseColor =
        eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    auto roughness =
        eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    auto anisotropic =
        eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    auto h = normalize(dir_in + dir_out);
    Spectrum Fresnel =
        baseColor + (1.0 - baseColor) * pow((1.0 - abs(dot(h, dir_out))), 5.0);
    Real aspect = sqrt(1.0 - anisotropic * 0.9);
    Real ax = max(0.001, roughness * roughness / aspect);
    Real ay = max(0.001, roughness * roughness * aspect);
    Real G = GTR(dir_in, ax, ay);
    Frame localFrame = Frame(h);
    Real hlx = dot(frame.x, h);
    Real hly = dot(frame.y, h);
    Real hlz = dot(frame.n, h);
    Real hlxyz = ((hlx * hlx) / pow(ax, 2.0) + (hlx * hlx) / pow(ax, 2.0) +
                  (hlx * hlx) / pow(ax, 2.0)) *
                 ((hlx * hlx) / pow(ax, 2.0) + (hlx * hlx) / pow(ax, 2.0) +
                  (hlx * hlx) / pow(ax, 2.0));
    Real D = 1.0 / (c_PI * ax * ay * hlxyz);
    Real pdf = D * G / (4.0 * abs(dot(vertex.shading_frame.n, dir_in)));
    return pdf;
    // Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Vector3 half_vector = normalize(dir_in + dir_out);
    // Real n_dot_in = dot(frame.n, dir_in);
    // Real n_dot_out = dot(frame.n, dir_out);
    // Real n_dot_h = dot(frame.n, half_vector);
    //
    // Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness);
    // Real D = GTR2(n_dot_h, roughness);
    // Real pdf = (G * D) / (4 * n_dot_in);
    // return pdf;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    auto baseColor =
        eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    auto roughness =
        eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    auto anisotropic =
        eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real aspect = sqrt(1.0 - anisotropic * 0.9);
    Real ax = max(0.001, roughness * roughness / aspect);
    Real ay = max(0.001, roughness * roughness * aspect);
    Spectrum Vh = normalize(Spectrum(ax * dir_in.x, ay * dir_in.y, dir_in.z));
    float lensq = Vh.x * Vh.x + Vh.y * Vh.y;
    Spectrum T1 = lensq > 0.0 ? Spectrum(-Vh.y, Vh.x, 0.0) * (1.0 / sqrt(lensq))
                              : Spectrum(1.0, 0.0, 0.0);
    Spectrum T2 = cross(Vh, T1);
    Real r = sqrt(rnd_param_uv.x);
    Real phi = 2.0 * c_PI * rnd_param_uv.y;
    Real t1 = r * cos(phi);
    Real t2 = r * sin(phi);
    Real s = 0.5 * (1.0 - Vh.z);
    t2 = (1.0 - s) * sqrt(1.0 - t1 * t1) + s * t2;
    Spectrum Nh =
        t1 * T1 + t2 * T2 + sqrt(max(0.0, 1.0 - t1 * t1 - t2 * t2)) * Vh;
    Spectrum Ne = normalize(Spectrum(ax * Nh.x, ay * Nh.y, max(0.0, Nh.z)));
    return BSDFSampleRecord{Ne, 0, roughness};

    // Vector3 local_dir_in = to_local(frame, dir_in);
    // Real roughness = eval(
    //     bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // roughness = std::clamp(roughness, Real(0.01), Real(1));
    // Real alpha = roughness * roughness;
    // Vector3 local_micro_normal =
    //     sample_visible_normals(local_dir_in, alpha, rnd_param_uv);
    //
    // Vector3 half_vector = to_world(frame, local_micro_normal);
    // Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    // return BSDFSampleRecord{
    //     reflected,
    //     Real(0) /* eta */, roughness /* roughness */
    // };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
