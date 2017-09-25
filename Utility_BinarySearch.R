#######
# This code is added by Tianchen Qian.
# The function binary.search() uses binary search to find 
# the root for a monotonic function.
#######

#######
# Major update log:
#     2017.08.01: implemented binary.search()
#     2017.09.06: implemented golden.section.search()
#######


# This function uses binary search to find the root for a monotonic function f.
# Input:
#     f: a monotonic function.
#     interval: starting interval of the binary search.
#     errtol: tolerance error.
# Output:
#     a number, which is the root of f.

binary.search <- function(f, interval = c(-20, 20), errtol = 1e-6) {
    if (length(interval) != 2) {
        stop("Interval in binarySearch should be vector of length 2 (lower, upper).")
    }
    lower <- interval[1]
    upper <- interval[2]
    
    if (f(lower) < f(upper)) {
        # f is an increasing function
        if ((f(lower) > 0) | (f(upper) < 0)) {
            stop("f is an increasing function, but f(lower) > 0 or f(upper) < 0.")
        }
        while(upper - lower > errtol) {
            mid <- (upper + lower) / 2
            if (f(mid) < 0) {
                lower <- mid
            } else {
                upper <- mid
            }
        }
    } else {
        # f is a decreasing function
        if ((f(lower) < 0) | (f(upper) > 0)) {
            stop("f is a decreasing function, but f(lower) < 0 or f(upper) > 0.")
        }
        while(upper - lower > errtol) {
            mid <- (upper + lower) / 2
            if (f(mid) < 0) {
                upper <- mid
            } else {
                lower <- mid
            }
        }
    }
    return(mid)
}




# This function uses binary search to minimize a concave function f.
# Input:
#     f: a concave function.
#     interval: starting interval of the binary search.
#     errtol: tolerance error.
# Output:
#     a number, which is the minimizer of f.

##### Implementing the golden section search method
##### a modification of the bisection method with the golden ratio
##### By Eric Cai - The Chemical Statistician
# copied from https://chemicalstatistician.wordpress.com/2013/04/22/using-r-to-implement-the-golden-bisection-method/

golden.section.search = function(f, lower.bound, upper.bound, tolerance = 1e-6)
{
    golden.ratio = 2/(sqrt(5) + 1)
    
    ### Use the golden ratio to set the initial test points
    x1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
    x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
    
    ### Evaluate the function at the test points
    f1 = f(x1)
    f2 = f(x2)
    
    iteration = 0
    
    while (abs(upper.bound - lower.bound) > tolerance)
    {
        iteration = iteration + 1
        # cat('', '\n')
        # cat('Iteration #', iteration, '\n')
        # cat('f1 =', f1, '\n')
        # cat('f2 =', f2, '\n')
        
        if (f2 > f1)
            # then the minimum is to the left of x2
            # let x2 be the new upper bound
            # let x1 be the new upper test point
        {
            # cat('f2 > f1', '\n')
            ### Set the new upper bound
            upper.bound = x2
            # cat('New Upper Bound =', upper.bound, '\n')
            # cat('New Lower Bound =', lower.bound, '\n')
            ### Set the new upper test point
            ### Use the special result of the golden ratio
            x2 = x1
            # cat('New Upper Test Point = ', x2, '\n')
            f2 = f1
            
            ### Set the new lower test point
            x1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
            # cat('New Lower Test Point = ', x1, '\n')
            f1 = f(x1)
        } 
        else 
        {
            # cat('f2 < f1', '\n')
            # the minimum is to the right of x1
            # let x1 be the new lower bound
            # let x2 be the new lower test point
            
            ### Set the new lower bound
            lower.bound = x1
            # cat('New Upper Bound =', upper.bound, '\n')
            # cat('New Lower Bound =', lower.bound, '\n')
            
            ### Set the new lower test point
            x1 = x2
            # cat('New Lower Test Point = ', x1, '\n')
            
            f1 = f2
            
            ### Set the new upper test point
            x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
            # cat('New Upper Test Point = ', x2, '\n')
            f2 = f(x2)
        }
    }
    
    ### Use the mid-point of the final interval as the estimate of the optimzer
    # cat('', '\n')
    # cat('Final Lower Bound =', lower.bound, '\n')
    # cat('Final Upper Bound =', upper.bound, '\n')
    estimated.minimizer = (lower.bound + upper.bound)/2
    # cat('Estimated Minimizer =', estimated.minimizer, '\n')
    print(estimated.minimizer)
}


# example
if (0) {
    golden.section.search( function(x){x^2}, -10, 10 )
    golden.section.search( function(x){x^2}, -1, 10 )
    golden.section.search( function(x){x^2}, -10, 1 )
    golden.section.search( function(x){x^2}, -10, -1 )
    golden.section.search( function(x){x^2}, 1, 10 )
}

