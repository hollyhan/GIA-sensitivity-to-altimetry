function md = initialize_model(md)
    
    % Initialize the model
    md.geometry.thickness = md.masstransport.spcthickness(1:end-1, 1);
    md.geometry.thickness(isnan(md.geometry.thickness)) = 0;
    md.geometry.surface(isnan(md.geometry.surface)) = 0;
    md.geometry.base(isnan(md.geometry.base)) = 0;
    md.mask.ice_levelset(isnan(md.mask.ice_levelset)) = 0;
    md.mask.ocean_levelset(isnan(md.mask.ocean_levelset)) = 0;
    md.geometry.surface = md.geometry.base + md.geometry.thickness;

end