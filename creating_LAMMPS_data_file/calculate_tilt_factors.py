import math

def calculate_tilt_factors(alpha, beta, gamma, a, b, c):
    # Convert angles from degrees to radians
    alpha_rad = math.radians(alpha)
    beta_rad = math.radians(beta)
    gamma_rad = math.radians(gamma)
    
    # Calculate the tilt factors
    xy = b * math.cos(gamma_rad)
    xz = c * math.cos(beta_rad)
    yz = c * (math.cos(alpha_rad) - math.cos(beta_rad) * math.cos(gamma_rad)) / math.sin(gamma_rad)
    
    return xy, xz, yz

def main():
    # Get user input for the angles
    alpha = float(input("Enter the angle α (in degrees): "))
    beta = float(input("Enter the angle β (in degrees): "))
    gamma = float(input("Enter the angle γ (in degrees): "))
    
    # Get user input for the lattice constants
    a = float(input("Enter the lattice constant a: "))
    b = float(input("Enter the lattice constant b: "))
    c = float(input("Enter the lattice constant c: "))
    
    # Calculate the tilt factors
    xy, xz, yz = calculate_tilt_factors(alpha, beta, gamma, a, b, c)
    
    # Print the tilt factors in the desired format
    print(f"{xy:.5f} {xz:.5f} {yz:.5f} xy xz yz")

if __name__ == "__main__":
    main()
