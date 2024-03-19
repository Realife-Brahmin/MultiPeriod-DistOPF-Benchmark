function PV = generatePVProfile(numTimePeriods, start, peak, finish)
    %GENERATEPVPROFILE Generate PV power coefficients over a time horizon.
    %
    % Synopsis:
    %   PV = generatePVProfile(numTimePeriods, start, peak, finish)
    %
    % Description:
    %   The function models the solar PV power generation over a specified
    %   time horizon using a cosine-based bell-shaped curve. The generation
    %   starts with a value of 'start', peaks at the 'peak' value, and 
    %   concludes with a value of 'finish'.
    %
    % Input:
    %   - numTimePeriods: Number of time intervals to consider.
    %   - start: PV generation value at the beginning.
    %   - peak: PV generation value at its peak.
    %   - finish: PV generation value at the end.
    %
    % Output:
    %   - PV: An array of size [1 x numTimePeriods] containing PV generation coefficients.
    %
    % Example:
    %   PV = generatePVProfile(1000, 0.7, 1.4, 0.5);
    %   % This will generate a PV profile starting at 0.7, peaking at 1.4, and finishing at 0.5.

    % Generate a time vector
    t = linspace(-pi/2, pi/2, numTimePeriods);
    
    % Calculate the PV coefficients using the cosine function and scaling
    PV = cos(t);
    
    % Rescale the coefficients based on start, peak, and finish values
    PV = start + (peak - start) * (PV + 1) / 2;
    
    midPoint = ceil(numTimePeriods/2);
    PV(midPoint:end) = linspace(peak, finish, floor(numTimePeriods/2) + 1);
end
