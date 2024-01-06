using Printf

function number_to_padded_string(num::Int)
    return @sprintf("%03d", num)
end

# Example usage:
println(number_to_padded_string(100))  # Outputs: "100"
println(number_to_padded_string(5))    # Outputs: "005"
println(number_to_padded_string(10))   # Outputs: "010"
