#include "../microfacet.h"

Real ior_fresnel(const Real eta,  // refracted / reflected ior
                  const Real kh)   // cosine of angle between normal/half-vector and direction
{
    Real costheta = 1.0f - (1.0f - kh * kh) / (eta * eta);

    if(costheta <= 0.0f)
    {
        return 1.0f;
    }

    costheta = sqrt(costheta);  // refracted angle cosine

    const Real n1t1 = kh;
    const Real n1t2 = costheta;
    const Real n2t1 = kh * eta;
    const Real n2t2 = costheta * eta;
    const Real r_p  = (n1t2 - n2t1) / (n1t2 + n2t1);
    const Real r_o  = (n1t1 - n2t2) / (n1t1 + n2t2);

    const Real fres = 0.5f * (r_p * r_p + r_o * r_o);

    return std::clamp(fres, 0.0, 1.0);
}

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    auto basecolor = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    auto roughness =
        eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    auto anisotropic =eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 h;
    if(reflect)
    {
        h = normalize(dir_in+dir_out);
    }else
    {
        h = normalize(dir_in+dir_out*eta);
    }
    if(dot(h,frame.n)<0.0)
    {
        h = -h;
    }
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real G = GTR(dir_in, roughness, anisotropic, frame) * GTR(
        dir_out, roughness, anisotropic, frame);
    Real aspect = sqrt(1.0 - anisotropic * 0.9);
    Real ax = max(0.001, roughness * roughness / aspect);
    Real ay = max(0.001, roughness * roughness * aspect);
    Real hlx = dot(frame.x, h);
    Real hly = dot(frame.y, h);
    Real hlz = dot(frame.n, h);
    Real hlxyz = ((hlx * hlx) / pow(ax, 2.0) + (hly * hly) / pow(ay, 2.0) +
                  (hlz * hlz)) *
                      ((hlx * hlx) / pow(ax, 2.0) + (hly * hly) / pow(ay, 2.0) +
                       (hlz * hlz));
    Real D = 1.0 / (c_PI * ax * ay * hlxyz);
    Real RS = (dot(h,dir_in)-eta*dot(h,dir_out))/(dot(h,dir_in)+eta*dot(h,dir_out));
    Real RP = (eta*dot(h,dir_in)-dot(h,dir_out))/(eta*dot(h,dir_in)+dot(h,dir_out));
    Real Fg = 0.5*(RS*RS+RP*RP);
    Fg = fresnel_dielectric(dot(h,dir_in), eta);
    if(reflect)
    {
        return (basecolor*G*D*Fg)/(4.0*abs(dot(frame.n,dir_in)));
    }
    return (sqrt(basecolor) * (1.0 - Fg) * D * G * abs(dot(h, dir_out) * dot(h, dir_in))) / (
        abs(dot(frame.n, dir_in))*pow((dot(h, dir_in) + eta * dot(h, dir_out)), 2.0));
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Vector3 h;
    if(reflect)
    {
        h = normalize(dir_in+dir_out);
    }else
    {
        h = normalize(dir_in+dir_out*eta);
    }
    if (dot(h, frame.n) < 0) {
        h = -h;
    }
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness,0.001,1.0);
    Real RS = (dot(h,dir_in)-eta*dot(h,dir_out))/(dot(h,dir_in)+eta*dot(h,dir_out));
    Real RP = (eta*dot(h,dir_in)-dot(h,dir_out))/(eta*dot(h,dir_in)+dot(h,dir_out));
    Real Fg = 0.5*(RS*RS+RP*RP);
    Real h_dot_in = dot(h, dir_in);
    Fg = fresnel_dielectric(h_dot_in, eta);
    auto anisotropic =eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real G = GTR(dir_in, roughness, anisotropic, frame) * GTR(
        dir_out, roughness, anisotropic, frame);
    Real aspect = sqrt(1.0 - anisotropic * 0.9);
    Real ax = max(0.001, roughness * roughness / aspect);
    Real ay = max(0.001, roughness * roughness * aspect);
    Real hlx = dot(frame.x, h);
    Real hly = dot(frame.y, h);
    Real hlz = dot(frame.n, h);
    Real hlxyz = ((hlx * hlx) / pow(ax, 2.0) + (hly * hly) / pow(ay, 2.0) +
                  (hlz * hlz)) *
                      ((hlx * hlx) / pow(ax, 2.0) + (hly * hly) / pow(ay, 2.0) +
                       (hlz * hlz));
    Real D = 1.0 / (c_PI * ax * ay * hlxyz);
    if(reflect)
    {
        return (Fg*D*G)/(4.0*fabs(dot(frame.n,dir_in)));
    }
    Real h_dot_out = dot(h,dir_out);
    Real sqrt_denom = h_dot_in+eta*h_dot_out;
    Real dh_dout = eta*eta*h_dot_out/(sqrt_denom*sqrt_denom);
    return (1.0-Fg)*D*G* fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
}


std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
   // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // Sample a micro normal and transform it to world space -- this is our half-vector.
    Real alpha = roughness * roughness;
    Vector3 local_dir_in = to_local(frame, dir_in);
    Vector3 local_micro_normal =
        sample_visible_normals(local_dir_in, alpha, rnd_param_uv);

    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    // Now we need to decide whether to reflect or refract.
    // We do this using the Fresnel term.
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);

    if (rnd_param_w <= F) {
        // Reflection
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        // set eta to 0 since we are not transmitting
        return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
    } else
    {
        // Refraction
        // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
        // (note that our eta is eta2 / eta1, and l = -dir_in)
        Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            // Total internal reflection
            // This shouldn't really happen, as F will be 1 in this case.
            return {};
        }
        // flip half_vector if needed
        if (h_dot_in < 0) {
            half_vector = -half_vector;
        }
        Real h_dot_out= sqrt(h_dot_out_sq);
        Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
        return BSDFSampleRecord{refracted, eta, roughness};
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
