function C_D = cd(mach,lam_C_f,tran_C_f,rough_C_f,C_f_C_M)
%% SETUP OF THE ROCKET'S GEOMETRY
N = 1; % number of fins
A_T = 1; % planform area of a single fin
A_r = 1; % reference area
f_b = 1; % fineness ratio of the body
A_W_B = 1; % total wetted area of the body
S = 1; % fin's span
gamma_L = 1; % leading edge sweepback angle
r_L = 1; % radius of the rounded leading edge
c_r = 1; % fin's root chord
h_r = 1; % root trailing edge thickness of the fin
t_r = 1; % maximum thickness of the fin
gamma_C = 1; % midcord line sweepback angle
A_B_f = 1; % fin's base area
f_n = 1; % fineness ratio of the nose
A_B = 1; % base area
%% C_D CALCULATION
if mach^2 <= 0.1 %incompressible flow
    %% COEFFICIENT OF FRICTION
    C_f = lam_C_f + tran_C_f + rough_C_f; % adding the C_f directly as just
    % one of them will be different than zero
    %% FRICTION DRAG COEFFICIENT
    C_D_f_t = 2 * N * C_f * A_T * A_r; % tail
    C_D_f_b = (1 + (0.5 / f_b)) * C_f * A_W_B / A_r; % body
    C_D_f = C_D_f_b + C_D_f_t; % friction drag coefficient
    %% FINS' DRAG COEFFICIENTS
    delta_C_D = ((1 - mach^2)^-0.417) - 1;
    C_D_L_T = 2 * N * S * r_L * cos(gamma_L)^2 * delta_C_D / A_r; % leading
    % edge pressure drag
    C_f_B = 2 * C_f * c_r / h_r;
    C_D_B_T = 0.135 * N / (C_f_B^(1/3)); % trailing edge base drag
    AR = S^2 / A_T; % aspect ratio of the fin
    epsilon = AR * t_r^(1/3) / c_r;
    C_D_TT_M = 1.15 * (t_r / c_r)^(5/3) * (1.61 + epsilon - ...
        sqrt((epsilon - 1.43)^2 + 0.578));
    K = cos(gamma_C)^2 + ((C_D_TT_M * A_r / (C_f_C_M * A_T) - ...
        4 * t_r * cos(gamma_C) / c_r) / ...
        (120 * (t_r / c_r)^4 * cos(gamma_C)^2));% prandtl correcting factor
    C_D_TT = 4 * N * C_f * (A_T / A_r) * (t_r * cos(gamma_C) / c_r + ...
        30 * (t_r / c_r)^4 * cos(gamma_C)^2 / ...
        ((K - mach^2 * cos(gamma_C)^2)^(3/2)));% thickness pressure drag
    % coefficient
    %% BODY'S DRAG COEFFICIENTS
    K_og = 1 + (6.82 * A_W_B * C_f_C_M * (f_n + 1)^2.22 / (f_b^3 * A_r));
    % ogive
    C_D_P = 6 * A_W_B * C_f / (f_b^3 * A_r * (K_og - mach^2)^0.6);% body
    % pressure drag
    C_D_B = 0.029 / sqrt(C_f * A_W_B / A_B); % body base drag
    %% TOTAL DRAG COEFFICIENT
    C_D = C_D_f + C_D_L_T + C_D_B_T + C_D_TT + C_D_P + C_D_B;
elseif mach <= 1 %subsonic flow
    %% CORRECTED COEFFICIENT OF FRICTION
    if lam_C_f ~= 0
        C_f_C = lam_C_f;
    elseif tran_C_f ~= 0
        C_f_C = tran_C_f * (1 - 0.09 * mach^2);
    else
        C_f_C = rough_C_f * (1 - 0.12 * mach^2);
    end
    %% FRICTION DRAG COEFFICIENT
    C_D_f_t = 2 * N * C_f_C * A_T * A_r; % tail
    C_D_f_b = (1 + (0.5/f_b)) * C_f_C * A_W_B / A_r; % body
    C_D_f = C_D_f_b + C_D_f_t; % friction drag coefficient
    %% FINS' DRAG COEFFICIENTS
    if mach < 0.9
        delta_C_D = ((1 - mach^2)^-0.417) - 1;
    else
        delta_C_D = 1 - 1.5 * (mach - 0.9);
    end
    C_D_L_T = 2 * N * S * r_L * cos(gamma_L)^2 * delta_C_D / A_r; % leading
    % edge pressure drag
    C_f_B = 2 * C_f_C * c_r / h_r;
    K_1 = cos(gamma_C)^2 + (0.223 + 4.02 * C_f_C_M * (t_r / h_r)^2)^2 / ...
        (c_r * C_f_C_M * h_r)^(2/3); % prandtl correcting factor
    C_D_B_T = 0.135 * N * A_B_f / (A_r * C_f_B^(1/3) * ...
        sqrt(K_1 - mach^2 * cos(gamma_C)^2)); % trailing edge base drag
    AR = S^2 / A_T; % aspect ratio of the fin
    epsilon = AR * t_r^(1/3) / c_r;
    C_D_TT_M = 1.15 * (t_r / c_r)^(5/3) * (1.61 + epsilon - ...
        sqrt((epsilon - 1.43)^2 + 0.578));
    K_2 = cos(gamma_C)^2 + ((C_D_TT_M * A_r / (C_f_C_M * A_T) - ...
        4 * t_r * cos(gamma_C) / c_r) / ...
        (120 * (t_r / c_r)^4 * cos(gamma_C)^2));% prandtl correcting factor
    C_D_TT = 4 * N * C_f_C * (A_T / A_r) * (t_r * cos(gamma_C) / c_r + ...
        30 * (t_r / c_r)^4 * cos(gamma_C)^2 / ...
        ((K_2 - mach^2 * cos(gamma_C)^2)^(3/2)));% thickness pressure drag
    % coefficient
    %% BODY'S DRAG COEFFICIENTS
    K_og = 1 + (6.82 * A_W_B * C_f_C_M * (f_n + 1)^2.22 / (f_b^3 * A_r));
    % ogive
    C_D_P = 6 * A_W_B * C_f_C / (f_b^3 * A_r * (K_og - mach^2)^0.6);% body
    % pressure drag
    K_3 = 1 + A_B / ((6.38 + 39.7 * h_r / c_r) * C_f_C_M * A_W_B);% prandtl
    % correcting factor
    C_D_B = 0.029 * A_B / ...
        (A_r * sqrt(C_f_C * A_W_B * (K_3 - mach^2) / A_B));
    %% TOTAL DRAG COEFFICIENT
    C_D = C_D_f + C_D_L_T + C_D_B_T + C_D_TT + C_D_P + C_D_B;
else %supersonic flow
    %% CORRECTED COEFFICIENT OF FRICTION
    if lam_C_f ~= 0
        C_f_C = lam_C_f / ((1 + 0.045 *mach^2)^0.25);
    elseif tran_C_f ~= 0
        C_f_C = tran_C_f / ((1 + 0.15 *mach^2)^0.58); % with K = 0.15
    else
        C_f_C = rough_C_f / (1 + 0.18 *mach^2);
    end
    %% FRICTION DRAG COEFFICIENT
    C_D_f_t = 2 * N * C_f_C * A_T * A_r; % tail
    C_D_f_b = (1 + (0.5 / f_b)) * C_f_C * A_W_B / A_r; % body
    C_D_f = C_D_f_b + C_D_f_t; % friction drag coefficient
    %% FINS' DRAG COEFFICIENTS
    delta_C_D = 1.214 - (0.502 / (mach^2)) + (0.1095 / (mach^4)) + ...
        (0.0231 / (mach^6));
    C_D_L_T = 2 * N * S * r_L * cos(gamma_L)^2 * delta_C_D / A_r; % leading
    % edge pressure drag
    C_D_B_T = N * (1 - 0.52 * mach^-1.19) * A_B_f / ...
        (A_r * (1 + 18 * C_f_C_M * (t_r / h_r)^2) * mach^2); % trailing
    % edge base drag
    %
    %% ADD HERE WAVE DRAG CALCULATION FOR THE FINS
    %
    %% BODY'S DRAG COEFFICIENTS
    %
    %% ADD HERE BODY PRESSURE DRAG CALCULATION
    %
    C_D_B_M = A_B * (0.185 + 1.15 * (h_r / c_r)) / A_r;
    M_cr = 0.892 * sqrt(C_D_B_M);
    if mach <= M_cr
        C_D_B = C_D_B_M * (0.88 + 0.12 * exp(-3.58 * (mach - 1)));
    else
        C_D_B = 0.7 * A_B / (A_r * mach^2);
    end
    %% TOTAL DRAG COEFFICIENT
    C_D = C_D_f + C_D_L_T + C_D_B_T + C_D_B; % wave drag and body pressure
    % drag contributions are missing
end