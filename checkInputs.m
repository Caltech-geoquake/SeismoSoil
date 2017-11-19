function [flag,err_msg] = checkInputs(data, option)


err_msg = string('');

if strcmpi(option(1:3),'all')
    [flag,err_msg] = checkAllInputs(data,option);
elseif strcmpi(option,'vs_profile')
    if ~isnumeric(data)
        flag = -1;
        err_msg = 'Data is not a numeric array. Check input file format.';
    elseif hasInfNaN(data)
        flag = -1;
        err_msg = 'Data has NaN.';
    elseif hasNegative(data)
        flag = -1;
        err_msg = 'Data has negative values.';
    elseif size(data,2) ~= 5  % number of columns
        flag = -1;
        err_msg = 'Data has incorrect number of columns. (Should be 5.)';
    elseif hasNonPositive(data(1:end-1,1))  % thickness
        flag = -1;
        err_msg = 'Soil thickness has non-positive numbers between first and second-to-last layers.';
    elseif hasNonPositive(data(:,2:4))  % Vs, damping, mass density
        flag = -1;
        err_msg = 'Vs, damping, or mass density has non-positive numbers.';
    elseif hasNonPositive(data(1:end-1,5))  % material number
        flag = -1;
        err_msg = 'Material number has non-positive numbers between first and second-to-last layers.';
    elseif any(data(:,3) > 1)  % damping larger than 1
        flag = -1;
        err_msg = 'Unit of damping ratio wrong: should be 1, not percent.';
    else
        flag = 1;
    end
elseif strcmpi(option,'curve')
    if ~isnumeric(data)
        flag = -1;
        err_msg = 'Data is not a numeric array. Check input file format.';
    elseif hasInfNaN(data)
        flag = -1;
        err_msg = 'Data has NaN.';
    elseif hasNegative(data)
        flag = -1;
        err_msg = 'Data has negative values.';
    elseif mod(size(data,2),4) ~= 0  % number of columns
        flag = -1;
        err_msg = 'Data has incorrect number of columns. (Should be a multiple of 4.)';
    elseif any(data(1:2:end,3) > 1)  % shear strain larger than 1
        flag = -1;
        err_msg = 'Unit of shear strain wrong: should be 1, not percent.';
    else
        flag = 1;
    end
elseif strcmpi(option,'motion')
    if ~isnumeric(data)
        flag = -1;
        err_msg = 'Data is not a numeric array. Check input file format.';
    elseif hasInfNaN(data)
        flag = -1;
        err_msg = 'Data has NaN.';
    elseif hasNegative(data(:,1))
        flag = -1;
        err_msg = 'Data has negative time values.';
    elseif size(data,2) ~= 2  % number of columns
        flag = -1;
        err_msg = 'Data has incorrect number of columns. (Should be 2.)';
    elseif ~ismonotonic(data(:,1),1,'INCREASING')
        flag = -1;
        err_msg = 'Data has time column not monotonically increasing.';
    else
        flag = 1;
    end
elseif any(strcmpi(option,{'h2n','h4g','h4x'}))
    if ~isnumeric(data)
        flag = -1;
        err_msg = 'Data is not a numeric array. Check input file format.';
    elseif hasInfNaN(data)
        flag = -1;
        err_msg = 'Data has NaN.';
    elseif hasNegative(data)
        flag = -1;
        err_msg = 'Data has negative values.';
    elseif size(data,1) ~= 4  % number of rows
        flag = -1;
        err_msg = 'Data has incorrect number of rows. (Should be 4.)';
    elseif any(data(2,:) ~= 0)
        flag = -1;
        err_msg = 'The second row should be all zeros.';
    else
        flag = 1;
    end
elseif any(strcmpi(option,{'hhg','hhx','hh_g','hh_x'}))
    if ~isnumeric(data)
        flag = -1;
        err_msg = 'Data is not a numeric array. Check input file format.';
    elseif hasInfNaN(data)
        flag = -1;
        err_msg = 'Data has NaN.';
    elseif hasNegative(data)
        flag = -1;
        err_msg = 'Data has negative values.';
    elseif size(data,1) ~= 9  % number of rows
        flag = -1;
        err_msg = 'Data has incorrect number of rows. (Should be 9.)';
    else
        flag = 1;
    end
elseif strcmpi(option,'tau_max')
    if ~isnumeric(data)
        flag = -1;
        err_msg = 'Data is not a numeric array. Check input file format.';
    elseif hasInfNaN(data)
        flag = -1;
        err_msg = 'Data has NaN.';
    elseif hasNegative(data)
        flag = -1;
        err_msg = 'Data has negative values.';
    elseif size(data,2) ~= 1  % number of columns
        flag = -1;
        err_msg = 'Data has incorrect number of rows. (Should be 1.)';
    else
        flag = 1;
    end
else
    error('Invalid "option" name.')
end

end

function [flag,err_msg] = checkAllInputs(data, option)

    err_msg = string('');
    
    if any(strcmpi(option,{'all_hh','all_h4'}))
        if strcmpi(option,'all_hh')
            id = 'HH';
        else
            id = 'H4';
        end
        vs_profile = data{1};
        curve = data{2};
        HHG = data{3};
        HHx = data{4};
        nr_materials = max(vs_profile(:,5));
        if size(HHG,2) ~= size(HHx,2)
            flag = -1;
            err_msg = sprintf('"%s_G" and "%s_x" dimensions do not match.',id,id);
        elseif size(curve,2) < nr_materials * 4
            flag = -1;
            err_msg = '"Curve" file does not have enough material properties.';
        elseif size(curve,2) ~= size(HHG,2)*4
            flag = -1;
            err_msg = sprintf('Number of columns in "curve" be 4 times of "%s_G" or "%s_x".',id,id);
        end
    elseif strcmpi(option,'all_h2')
        vs_profile = data{1};
        curve = data{2};
        H2n = data{3};
        nr_materials = max(vs_profile(:,5));
        if size(curve,2) < nr_materials * 4
            flag = -1;
            err_msg = '"Curve" file does not have enough material properties.';
        elseif size(curve,2) ~= size(H2n,2)*4
            flag = -1;
            err_msg = 'Number of columns in "curve" be 4 times of "H2n".';
        end
    elseif strcmpi(option,'all_eql')
        vs_profile = data{1};
        curve = data{2};
        nr_materials = max(vs_profile(:,5));
        if size(curve,2) < nr_materials * 4
            flag = -1;
            err_msg = '"Curve" file does not have enough material properties.';
        end
    else
        error('Invalid "option" name.');
    end
end


function tf = hasNegative(x)
    tf = any(x(:) < 0);
end

function tf = hasNonPositive(x)
    tf = any(x(:) <= 0);
end

