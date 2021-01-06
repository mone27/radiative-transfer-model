## Porting of the MatLab code from 14_02

# Supplemental program 14.2

# ---------------------------------------------------------
  # Test two-stream model with clumping. Compare numerical
# integration of leaf fluxes with canopy integrated fluxes.
# ---------------------------------------------------------
  
# Input parameters

  
# --- Leaf optical properties

rho_leaf  <- 0.10;                   # Leaf reflectance
tau_leaf <- 0.05;                    # Leaf transmittance
omega_leaf <- rho_leaf + tau_leaf;   # Leaf scattering coefficient
Kb <- 0.58;                          # Direct beam extinction coefficient
Kd <- 0.70;                          # Diffuse extinction coefficient
beta <- 0.54;                        # Upscatter parameter for diffuse radiation
beta0 <- 0.46;                       # Upscatter parameter for direct beam radiation

# --- Canopy variables

LAI <- 6;                            # Leaf area index (m2/m2)
OMEGA <- 0.75;                       # Clumping index

# --- Soil variables

albsoib <- 0.1;                      # Soil albedo (direct beam)
albsoid <- 0.1;                      # Soil albedo (diffuse)

# --- Atmospheric solar radiation, given here as a unit of incoming radiation

swskyb <- 0.8;                       # Direct beam (W/m2)
swskyd <- 0.2;                       # Diffuse (W/m2)


# Analytical solution

  
# --- Common terms: Eqs. (14.87) - (14.91)

b <- (1 - (1 - beta) * omega_leaf) * Kd;
c <- beta * omega_leaf * Kd;
h <- sqrt(b*b - c*c);
u <- (h - b - c) / (2 * h);
v <- (h + b + c) / (2 * h);
g1 <- (beta0 * Kb - b * beta0 - c * (1 - beta0)) * omega_leaf * Kb * swskyb / (h*h - Kb*Kb);
g2 <- ((1 - beta0) * Kb + c * beta0 + b * (1 - beta0)) * omega_leaf * Kb * swskyb / (h*h - Kb*Kb);

# --- Exponential functions of leaf area

s1 <- function(x) exp(-h * OMEGA * x);
s2 <- function(x) exp(-Kb * OMEGA * x);

# --- Direct beam solution

# n1 (Eq. 14.92) and n2 (14.93)

num1 <- v * (g1 + g2 * albsoid + albsoib * swskyb) * s2(LAI);
num2 <- g2 * (u + v * albsoid) * s1(LAI);
den1 <- v * (v + u * albsoid) / s1(LAI);
den2 <- u * (u + v * albsoid) * s1(LAI);
n2b <- (num1 - num2) / (den1 - den2);
n1b <- (g2 - n2b * u) / v;

# Scattered direct beam fluxes:
# iupwb - direct beam flux scattered upward above cumulative LAI (W/m2); Eq. (14.94)
# idwnb - direct beam flux scattered downward below cumulative LAI (W/m2); Eq. (14.95)
# and their derivatives with respect to LAI

iupwb <-  function(x) -g1 * s2(x) + n1b * u * s1(x) + n2b * v / s1(x);
idwnb <-  function(x)  g2 * s2(x) - n1b * v * s1(x) - n2b * u / s1(x);
diupwb <- function(x) ( Kb * g1 * s2(x) - h * n1b * u * s1(x) + h * n2b * v / s1(x)) * OMEGA;
didwnb <- function(x) (-Kb * g2 * s2(x) + h * n1b * v * s1(x) - h * n2b * u / s1(x)) * OMEGA;

# icb - direct beam flux absorbed by canopy (W/m2); Eq. (14.97)

icb <- swskyb * (1 - s2(LAI)) - iupwb(0) + iupwb(LAI) - idwnb(LAI);

# icsunb - direct beam flux absorbed by sunlit canopy (W/m2); Eq. (14.114)
# icshab - direct beam flux absorbed by shaded canopy (W/m2); Eq. (14.115)

a1b <- -g1 *      (1 - s2(LAI)*s2(LAI)) / (2 * Kb) + 
n1b * u * (1 - s2(LAI)*s1(LAI)) / (Kb + h) + n2b * v * (1 - s2(LAI)/s1(LAI)) / (Kb - h);
a2b <-  g2 *      (1 - s2(LAI)*s2(LAI)) / (2 * Kb) -
n1b * v * (1 - s2(LAI)*s1(LAI)) / (Kb + h) - n2b * u * (1 - s2(LAI)/s1(LAI)) / (Kb - h);

icsunb <- (1 - omega_leaf) * ((1 - s2(LAI)) * swskyb + Kd * (a1b + a2b) * OMEGA);
icshab <- icb - icsunb;


# --- Diffuse solution

# n1 (Eq. 14.99) and n2 (14.100)

num <- swskyd * (u + v * albsoid) * s1(LAI);
den1 <- v * (v + u * albsoid) / s1(LAI);
den2 <- u * (u + v * albsoid) * s1(LAI);
n2d <- num / (den1 - den2);
n1d <- -(swskyd + n2d * u) / v;

# Scattered diffuse fluxes:
# iupwd - diffuse flux scattered upward above cumulative LAI (W/m2); Eq. (14.101)
# idwnd - diffuse flux scattered downward below cumulative LAI (W/m2); Eq. (14.102)
# and their derivatives with respect to LAI

iupwd <- function(x)  n1d * u * s1(x) + n2d * v / s1(x);
idwnd <- function(x) -n1d * v * s1(x) - n2d * u / s1(x);
diupwd <- function(x) (-h * n1d * u * s1(x) + h * n2d * v / s1(x)) * OMEGA;
didwnd <- function(x) ( h * n1d * v * s1(x) - h * n2d * u / s1(x)) * OMEGA;

# icd - diffuse flux absorbed by canopy (W/m2); Eq. (14.104)

icd <- swskyd - iupwd(0) + iupwd(LAI) - idwnd(LAI);

# icsund - diffuse flux absorbed by sunlit canopy (W/m2); Eq. (14.118)
# icshad - diffuse flux absorbed by shaded canopy (W/m2); Eq. (14.119)

a1d <-  n1d * u * (1 - s2(LAI)*s1(LAI)) / (Kb + h) + n2d * v * (1 - s2(LAI)/s1(LAI)) / (Kb - h);
a2d <- -n1d * v * (1 - s2(LAI)*s1(LAI)) / (Kb + h) - n2d * u * (1 - s2(LAI)/s1(LAI)) / (Kb - h);

icsund <- (1 - omega_leaf) * Kd * (a1d + a2d) * OMEGA;
icshad <- icd - icsund;

# --- Total canopy flux

ic <- icb + icd;
icsun <- icsunb + icsund;
icsha <- icshab + icshad;

# ---------------------------
  # Print output
# ---------------------

ic
icsun
icsha

