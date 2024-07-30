from raytracing.spherical_grin import raytrace_spherical
from raytracing.axial_grin import raytrace_axial
from raytracing.cylindrical_grin import raytrace_cylindrical

def main():
    # Example parameters for the raytracing functions
    params_spherical = {...}  # Replace with actual parameters
    params_linear = {...}  # Replace with actual parameters
    params_cylindrical = {...}  # Replace with actual parameters

    # Run the raytracing functions
    result_spherical = raytrace_spherical(params_spherical)
    result_linear = raytrace_axial(params_linear)
    result_cylindrical = raytrace_cylindrical(params_cylindrical)

    # Print or process the results
    print("Spherical Raytracing Result:", result_spherical)
    print("Linear Raytracing Result:", result_linear)
    print("Cylindrical Raytracing Result:", result_cylindrical)

if __name__ == "__main__":
    main()
