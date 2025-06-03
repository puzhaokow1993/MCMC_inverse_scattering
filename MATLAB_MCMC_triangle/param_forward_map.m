function far = param_forward_map(k, thetaout, theta, coefficient_F,XY_coordinates)
    % Ensure coefficient_F is structured correctly
    % coefficient_F should be a 3-by-N array:
    % - First row: x-coordinates
    % - Second row: y-coordinates
    % - Third row: function values F(x, y)

    % Extract components
    x = XY_coordinates(1, :);
    y = XY_coordinates(2, :);
    n = link_func(coefficient_F);

    % Create the scatteredInterpolant object once
    F = scatteredInterpolant(x', y', n', 'natural', 'none');

    % Create a function handle for evaluation
    interpFunc = @(xq, yq) F(xq, yq);

    % Call the forward_map function
    far = forward_map(k, thetaout, theta, interpFunc);
end
