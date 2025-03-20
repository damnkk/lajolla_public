inline auto Fss(const eval_op* eval_op,const DisneyDiffuse& bsdf,const Spectrum& v,const Frame& frame)-> Real
{
    Spectrum h = normalize(eval_op->dir_in +eval_op-> dir_out);
    Real Fss90 = eval(bsdf.roughness,eval_op->vertex.uv,eval_op->vertex.uv_screen_size,eval_op->texture_pool)*abs(dot(h,eval_op->dir_out)*dot(h,eval_op->dir_out));
    return (1.0+(Fss90-1.0)*pow(1.0-abs(dot(frame.n,v)),5));
}

inline auto FD(const eval_op* eval_op, const DisneyDiffuse& bsdf, const Spectrum& v,
               const Frame& frame) -> Real
{
    Spectrum h = normalize(eval_op->dir_in +eval_op-> dir_out);
    Real Fd90 = 0.5+2.0*eval(bsdf.roughness,eval_op->vertex.uv,eval_op->vertex.uv_screen_size,eval_op->texture_pool)*abs(dot(h,eval_op->dir_out)*dot(h,eval_op->dir_out));
    return (1.0+(Fd90-1.0)*pow(1.0-abs(dot(frame.n,v)),5));
}
Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    Spectrum diffuseBsdf =FD(this,bsdf,dir_in,frame)*FD(this,bsdf,dir_out,frame)* fmax(dot(frame.n, dir_out), Real(0)) *
          eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool) / c_PI;

    Real FssIn = Fss(this,bsdf,dir_in,frame);
    Real FssOut = Fss(this,bsdf,dir_out,frame);
    Real absNdotI = abs(dot(frame.n,dir_in));
    Real absNdotO = abs(dot(frame.n,dir_out));
    Spectrum diffuseSubsurfaceBsdf =(FssIn*FssOut*(1.0/(absNdotI+absNdotO)-0.5)+0.5)*absNdotO* 1.25*eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool) / c_PI;
    // Homework 1: implement this!
    return (1.0- eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool))*diffuseBsdf+eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool)*diffuseSubsurfaceBsdf;
    return diffuseBsdf;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(1) /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
