function result = generateLoadProfile(n, low, high, trend)
    % Given a number of time periods (n), the function generates an array,
    % which starts and ends at the provided low value, and goes to the high value in the middle.
    % Default low and high values are 1.0 and 1.4, respectively.

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
    if mod(n,2) == 0
        half_n = n/2;
        first_half = linspace(low, high, half_n);
        second_half = linspace(high, low, half_n);
        result = [first_half, second_half];
    else
        half_n = (n-1)/2;
        first_half = linspace(low, high, half_n + 1);
        second_half = linspace(high, low, half_n);
        result = [first_half, second_half(2:end)];
    end
    
end
