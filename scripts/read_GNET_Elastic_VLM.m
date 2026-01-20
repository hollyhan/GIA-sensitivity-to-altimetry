%% read_GNET_Elastic_VLM.m
% Reads Table S1 of Berg et al., 2023 (https://doi.org/10.1029/2023GL104851)
%
% Expected columns:
%   Station | tstart [yr] | tend [yr] |
%   Uobserved [mm/yr] | Uelastic,GrIS [mm/yr] | 
%   Uelastic,GrPG [mm/yr] | Uelastic,CanPG [mm/yr] | UGIA,obs [mm/yr]
%
% Usage:
%   data = read_GNET_Elastic_VLM('Table_S1_GNET_VLM_Berg_et_al.xlsx', true);

function data = read_GNET_VLM_table(filename, plot_fig)

    if nargin < 1
        filename = 'Table_S1_GNET_VLM_Berg_et_al.xlsx';
    end
    if nargin < 2
        plot_fig = true;
    end

    % --- Read the Excel table ---
    fprintf('\n📘 Reading %s ...\n', filename);
    opts = detectImportOptions(filename, 'VariableNamingRule','preserve');
    tbl = readtable(filename, opts);

    % --- Parse into structure ---
    data.station  = string(tbl.Station);
    data.tstart   = tbl.("tstart [yr]");
    data.tend     = tbl.("tend [yr]");

    % --- Handle the "value ± sigma" string fields ---
    numericFields = ["Uobserved [mm/yr]", ...
                     "Uelastic,GrIS [mm/yr]", ...
                     "Uelastic,GrPG [mm/yr]", ...
                     "Uelastic,CanPG [mm/yr]", ...
                     "UGIA,obs [mm/yr]"];

    for f = 1:numel(numericFields)
        col = tbl.(numericFields(f));
        val = nan(height(tbl),1);
        sig = nan(height(tbl),1);

        for i = 1:height(tbl)
            str = strrep(string(col{i}), ' ', '');
            parts = split(str, '±');
            if numel(parts) == 2
                val(i) = str2double(parts(1));
                sig(i) = str2double(parts(2));
            elseif numel(parts) == 1
                val(i) = str2double(parts(1));
                sig(i) = NaN;
            end
        end

        % --- Clean field names ---
        fieldName = char(matlab.lang.makeValidName(numericFields(f)));
        fieldName = regexprep(fieldName, '_mm_yr_?$', '', 'ignorecase');  % remove trailing underscores

        data.(fieldName) = val;
        data.([fieldName '_sigma']) = sig;
    end

    % --- Print summary ---
    fprintf('✅ Successfully read %d stations.\n', height(tbl));
    fprintf('  GNSS Time span varies between: %.1f – %.1f yr\n', min(data.tstart), max(data.tend));

    % --- Optional diagnostic plots ---
    if plot_fig
        figure('Color','w','Position',[400 300 950 500]);
        bar(mean([data.Uelastic_GrIS, ...
                  data.Uelastic_GrPG, ...
                  data.Uelastic_CanPG], 'omitnan'));
        set(gca, 'XTickLabel', {'GrIS','GrPG','CanPG'});
        ylabel('Average elastic uplift (mm/yr)');
        title('Mean Elastic Contributions');
    end
end
