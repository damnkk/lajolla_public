#include <complex>

#include "../microfacet.h"

Real R0(Real param)
{
    return pow(param-1.0,2.0)/pow(param+1.0,2.0);
}
Real GGX_clearCoat(const Spectrum& dir,const Frame& frame)
{
    Spectrum  wl= to_local(frame,dir);
    Real A = (sqrt(1.0+(pow(wl.x*0.25,2.0)+pow(wl.y*0.25,2.0))/pow(wl.z,2.0))-1.0)/2.0;
    return 1.0/(1.0+A);
}
Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    // Homework 1: implement this!
    Spectrum h= normalize(dir_in+dir_out);
    Real clearcoatGloss = eval(bsdf.clearcoat_gloss,vertex.uv,vertex.uv_screen_size,texture_pool);
    Real Fc = R0(1.5)+(1.0-R0(1.5))*pow(1.0-abs(dot(frame.n,dir_out)),5);
    Real ag = (1.0-clearcoatGloss)*0.1+clearcoatGloss*0.001;
    Real D = (ag*ag-1.0)/(c_PI*log(ag*ag)*(1.0+(ag*ag-1.0)*pow(to_local(frame,h).z,2.0)));
    Real Gc = GGX_clearCoat(dir_out,frame)*GGX_clearCoat(dir_in,frame);
    return make_const_spectrum( Fc*Gc*D/(4.0*abs(dot(frame.n,dir_in))));
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Spectrum h= normalize(dir_in+dir_out);
    Real clearcoatGloss = eval(bsdf.clearcoat_gloss,vertex.uv,vertex.uv_screen_size,texture_pool);
    Real ag = (1.0-clearcoatGloss)*0.1+clearcoatGloss*0.001;
    Real D = (ag*ag-1.0)/(c_PI*log(ag*ag)*(1.0+(ag*ag-1.0)*pow(to_local(frame,h).z,2.0)));
    Real pdfCoat = D*dot(frame.n,h)/(4.0*abs(dot(dir_out,h)));
    return pdfCoat;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoatGloss = eval(bsdf.clearcoat_gloss,vertex.uv,vertex.uv_screen_size,texture_pool);
    Real ag = (1.0-clearcoatGloss)*0.1+clearcoatGloss*0.001;
    Real cosH = sqrt((1.0-pow(pow(ag,2.0),1-rnd_param_uv.x))/(1.0-pow(ag,2.0)));
    Real h_elevation = acos(cosH);
    Real sinH = sin(h_elevation);
    Real Hazimuth = 2.0*c_PI*rnd_param_uv.y;
    Real hlx = sinH*cos(Hazimuth);
    Real hly = sinH*sin(Hazimuth);
    Real hlz = cosH;
    Vector3 half_vector = to_world(frame,Spectrum(hlx,hly,hlz));
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{reflected, 0, ag};
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}

