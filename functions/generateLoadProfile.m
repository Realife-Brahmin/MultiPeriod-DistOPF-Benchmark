function result = generateLoadProfile(n, low, high, trend)
    %GENERATELOADPROFILE Generates a symmetrical load profile.
    % This function generates a load profile over `n` time periods that starts
    % and ends at the provided low value and peaks at the high value in the middle.
    % The profile is symmetrical about its midpoint.
    %
    % Syntax:
    %   result = generateLoadProfile(n, low, high, trend)
    %
    % Inputs:
    %   n     - The number of time periods for the load profile.
    %   low   - Optional. The starting and ending value of the profile. Default is 1.0.
    %   high  - Optional. The peak value of the profile, occurring at the midpoint. Default is 1.4.
    %   trend - Optional. A string specifying the profile shape. Currently supports:
    %           'normal'  - Symmetrical peak profile (default).
    %           'upwards' - Linearly increasing profile from `low` to `high`.
    %
    % Outputs:
    %   result - A 1-by-n array containing the generated load profile.
    %
    % Examples:
    %   % Generate a symmetrical load profile for 10 periods:
    %   profile = generateLoadProfile(10);
    %
    %   % Generate a load profile with specific low and high values:
    %   profile = generateLoadProfile(10, 0.8, 1.5);
    %
    %   % Generate a linearly increasing load profile:
    %   profile = generateLoadProfile(10, 0.5, 1.5, 'upwards');

    if nargin < 2
        low = 1.0;
    end

    if nargin < 3
        high = 1.4;
    end

    if nargin < 4
        trend = 'normal';
    end
    
    if strcmp(trend, 'upwards')
        result = linspace(low, high, n);
        return;
    end

    % Check if n is even or odd
    if mod(n, 2) == 0
        % For even n, ensure that the two middle values are 'high'
        half_n = n / 2;
        first_half = linspace(low, high, half_n);
        second_half = linspace(high, low, half_n);
        result = [first_half, second_half]; % Ensures middle two values are 'high'
    else
        % For odd n, ensure the exact middle value is 'high'
        half_n = (n - 1) / 2;
        first_half = linspace(low, high, half_n + 1);
        second_half = linspace(high, low, half_n + 1);
        result = [first_half, second_half(2:end)]; % Removes the duplicated 'high' at the middle
    end
    
end
