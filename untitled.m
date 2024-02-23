findLargestDigit(12345)

function largestDigit = findLargestDigit(x)
% finds the largest digit in a given integer
%   INPUT: x -  is an integer
%  RETURN: largestDigit -  if x is an integer, digit between 0 and 9 in x
%                          if x is NOT an integer, returns -1

% initialize
largestDigit = 0;
if rem(x,1) ~= 0 % not an integer
    largestDigit = -1;
else % it is an integer
    while abs(x) > 0
        % pick the last digit
        digit = abs(rem(x, 10));
        
        % check if it is the largest
        if digit > largestDigit
            largestDigit = digit;
        end

        % update x to get rid of the last digit
        x = (x - rem(x, 10))/10;
    end
end
end