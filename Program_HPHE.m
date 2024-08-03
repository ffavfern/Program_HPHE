function number_tube = Program_HPHE(Eff_c_opti, Eff_thermosy)

    % แสดงข้อความเริ่มต้นโปรแกรม
    fprintf('เริ่มทำงาน Program_HPHE...\n');
    
    % ตรวจสอบและกำหนดค่าเริ่มต้นให้กับตัวแปรถ้าไม่มีการป้อนค่าเข้ามา
    if nargin < 2
        Eff_c_opti = 0.9;
        Eff_thermosy = 0.8;
    end
    
    % แสดงค่าเริ่มต้นของประสิทธิภาพ
    fprintf('Eff_c_opti: %.2f, Eff_thermosy: %.2f\n', Eff_c_opti, Eff_thermosy);

    % กำหนดจำนวนรอบการทำงานสูงสุดของลูป
    maxIterations = 1000; 
    iterationCount = 0;
    error_C_Eff = abs(Eff_c_opti - Eff_thermosy);

    % ฟังก์ชันย่อยสำหรับการสุ่มหมายเลขท่อ
    function newTubeNumber = getRandomTubeNumber()
        fprintf('กำลังสุ่มหมายเลขท่อ...\n');
        newTubeNumber = randi([1, 480]); 
        fprintf('หมายเลขท่อที่เลือก: %d\n', newTubeNumber);
    end

    % กำหนดค่าเริ่มต้นของจำนวนท่อ
    n_tube = getRandomTubeNumber();

    % สร้างอาเรย์สำหรับเก็บค่าประสิทธิภาพและค่าประหยัดสุทธิ
    efficiency_values = [];
    net_savings_values = [];
    
    % วนลูปจนกว่าค่าความผิดพลาดจะต่ำกว่าค่าที่กำหนด
    while error_C_Eff > 0.009
        iterationCount = iterationCount + 1;

        % หยุดลูปหากจำนวนรอบเกินค่าที่กำหนด
        if iterationCount > maxIterations
            error('ถึงจำนวนรอบสูงสุดแล้ว หยุดการทำงานเพื่อป้องกันการทำงานที่ไม่มีที่สิ้นสุด');
        end

        % แสดงค่าความผิดพลาดและรอบการทำงานในแต่ละรอบ
        fprintf('รอบที่: %d, error_C_Eff: %.6f\n', iterationCount, error_C_Eff);

        try
            fprintf('เรียกใช้ฟังก์ชัน Optimum_HPHE ด้วย n_tube: %d\n', n_tube);
            [Eff_thermosy, Eff_c_opti, SL, column, n_row2, ST, Le, Lc, P1, Q1, P2, CA, Cc, C_st, U, CE, H, m_cp_min, T_hi, T_ci] = Optimum_HPHE(n_tube);

            % ตรวจสอบค่า column และ n_row2 ว่าเป็นค่าที่ถูกต้อง
            if column <= 0 || n_row2 <= 0
                error('ค่าที่ไม่ถูกต้องสำหรับ column หรือ n_row2: column = %d, n_row2 = %d', column, n_row2);
            end

            % คำนวณค่า n_row และปรับค่า n_tube2 ตามค่า n_row
            n_row = min(max(floor(n_tube / column), 1), 9);  
            n_tube2 = n_row * column;  

            fprintf('คำนวณค่า n_row เป็น %d\n', n_row);

            % คำนวณค่าความยาวเชิงเส้น (SD)
            SD = sqrt(ST^2 + SL^2);
            fprintf('คำนวณค่า SD เป็น %.4f\n', SD);

            % ตรวจสอบค่า SD ว่าเป็นค่าที่ถูกต้อง
            if isnan(SD) || isinf(SD)
                error('ค่าที่ไม่ถูกต้องสำหรับ SD: %f', SD);
            end

            % คำนวณค่าประสิทธิภาพและค่าประหยัดสุทธิ
            ii = 1;
            EFF_c1 = [];
            S = [];
            for eff = 0.01:0.02:0.98
                EFF_c1(ii) = eff;
                non = 1 - EFF_c1(ii) * C_st;
                if non >= 0
                    S(ii) = P1 * CE * H * EFF_c1(ii) * m_cp_min * (T_hi - T_ci) - P2 * CA * Cc * ...
                        log((1 - EFF_c1(ii) * C_st) / (1 - EFF_c1(ii))) / ((1 - C_st) * U);
                    efficiency_values(ii) = EFF_c1(ii);
                    net_savings_values(ii) = S(ii);
                end
                ii = ii + 1;
            end

            % อัปเดตค่าความผิดพลาดและแสดงผลค่าที่เปลี่ยนแปลง
            previous_error_C_Eff = error_C_Eff;
            error_C_Eff = abs(Eff_c_opti - Eff_thermosy);
            fprintf('ค่าความผิดพลาดก่อนหน้า: %.6f, ค่าความผิดพลาดหลังจากอัปเดต: %.6f\n', previous_error_C_Eff, error_C_Eff);

            % หยุดลูปหากค่าความผิดพลาดเปลี่ยนแปลงไม่มาก
            if abs(previous_error_C_Eff - error_C_Eff) < 1e-6
                fprintf('ค่าความผิดพลาดไม่เปลี่ยนแปลงมาก หยุดลูปเพื่อป้องกันการทำงานที่ไม่มีที่สิ้นสุด\n');
                break;
            end

            % สุ่มหมายเลขท่อใหม่สำหรับรอบถัดไป
            n_tube = getRandomTubeNumber();  

        catch ME
            % แสดงข้อผิดพลาดที่เกิดขึ้น
            fprintf('เกิดข้อผิดพลาด: %s\n', ME.message);
            if strcmp(ME.identifier, 'MATLAB:dimagree')
                error('ขนาดของเมทริกซ์ต้องเท่ากัน: %s', ME.message);
            else
                rethrow(ME);
            end
        end
    end

    % คืนค่าหมายเลขท่อที่ได้จากการคำนวณ
    number_tube = n_tube2;  
    fprintf('Program_HPHE เสร็จสิ้นการทำงาน หมายเลขท่อสุดท้ายคือ: %d\n', number_tube);
    
    % --- แสดงผลลัพธ์ ---
    fprintf('\n--- ผลลัพธ์สุดท้าย ---\n');
    fprintf('ความยาวของ Evaporator (Le): %.2f m\n', Le);
    fprintf('ความยาวของ Condenser (Lc): %.2f m\n', Lc);
    fprintf('จำนวนท่อแลกเปลี่ยนความร้อน: %d\n', number_tube);
    fprintf('ประสิทธิภาพทางเศรษฐศาสตร์ที่เหมาะสมที่สุด: %.4f\n', Eff_c_opti);

    % คำนวณระยะเวลาคืนทุน (สมมติว่า P1 คือระยะเวลาคืนทุน)
    payback_period = P1; 
    fprintf('ระยะเวลาคืนทุน: %.2f ปี\n', payback_period);

    % คำนวณมูลค่าประหยัดสุทธิ (สมมติว่า Q1 เป็นค่าพลังงานที่ประหยัดได้)
    net_savings = Q1; 
    fprintf('มูลค่าประหยัดสุทธิ: %.2f หน่วยเงิน\n', net_savings);

    % แสดงค่าประสิทธิภาพทางเศรษฐศาสตร์ที่ดีที่สุด
    best_eff_economic = Eff_c_opti; 
    fprintf('ประสิทธิภาพทางเศรษฐศาสตร์ที่ดีที่สุด: %.4f\n', best_eff_economic);

    % สร้างกราฟแสดงผลหลังจากการคำนวณเสร็จสิ้น
    if ~isempty(efficiency_values) && ~isempty(net_savings_values)
        figure;
        plot(efficiency_values, net_savings_values, 'r-', 'DisplayName', 'Net Savings');
        hold on;
        plot(efficiency_values, efficiency_values, 'b--', 'DisplayName', 'Efficiency Values');
        hold off;
        axis tight;
        title('Effectiveness vs Value');
        xlabel('Effectiveness');
        ylabel('Value');
        legend('show');
        fprintf('แสดงผลลัพธ์และกราฟเรียบร้อยแล้ว\n');
    else
        fprintf('ไม่มีข้อมูลสำหรับสร้างกราฟ\n');
    end
end

function [Eff_thermosy, Eff_c_opti, Q1, P1, CE, H, m_cp_min, T_hi, T_ci, P2, CA, Cc, C_st, U, i, d, A_HX, N_day, column, n_row2] = Optimum_HPHE(n_tube)
    try
        % แสดงข้อความเริ่มต้นการทำงานของฟังก์ชัน
        fprintf('เริ่มต้นฟังก์ชัน Optimum_HPHE ด้วย n_tube = %d\n', n_tube);

        % รับอินพุตที่เป็นค่าคงที่จากผู้ใช้หรือการคำนวณเริ่มต้น
        Le = 0.345; % ความยาวของส่วน Evaporator (เมตร)
        Lc = 0.345; % ความยาวของส่วน Condenser (เมตร)
        di = 0.01905; % เส้นผ่านศูนย์กลางภายในของท่อ (เมตร)
        do = 0.01921; % เส้นผ่านศูนย์กลางภายนอกของท่อ (เมตร)
        ST = 0.04; % ระยะพิทช์ตามขวาง (เมตร)
        SL = 0.04; % ระยะพิทช์ตามยาว (เมตร)
        T_hi = 200 + 273.15; % อุณหภูมิของของไหลที่เข้าให้ความร้อน (เคลวิน)
        T_ci = 52 + 273.15; % อุณหภูมิของของไหลที่เข้าระบายความร้อน (เคลวิน)
        Q_h = 576; % อัตราการไหลของของไหลที่ให้ความร้อน (ลูกบาศก์เมตร/ชั่วโมง)
        Q_co = 0.33; % อัตราการไหลของของไหลที่ระบายความร้อน (ลูกบาศก์เมตร/ชั่วโมง)
        column = 10; % จำนวนคอลัมน์ของท่อ
        heating_medium = 2; % ตัวกลางให้ความร้อน; 1 สำหรับน้ำ, 2 สำหรับอากาศ
        cooling_medium = 1; % ตัวกลางระบายความร้อน; 1 สำหรับน้ำ, 2 สำหรับอากาศ
        characteristic_Of_HPHE = 1; % การจัดเรียงท่อ; 1 สำหรับ Staggered, 2 สำหรับ Aligned
        materials = 2; % วัสดุที่ใช้ทำท่อ; 1 สำหรับเหล็ก, 2 สำหรับทองแดง, 3 สำหรับอลูมิเนียม
        Working_fluid = 1; % ของไหลที่ใช้ในการทำงาน; 1 สำหรับน้ำ, 2 สำหรับ R134a
        fillratio = 0.5; % อัตราการเติมของของไหลที่ใช้ในการทำงาน
        bata = 1.57; % มุมเอียง
        Df_e = 0; % ส่วนของ Evaporator
        Df_c = 0; % ส่วนของ Condenser

        % คำนวณค่าพื้นฐานของการแลกเปลี่ยนความร้อน
        n_tube2 = max(n_tube, column);
        n_tube2 = ceil(n_tube2 / column) * column; % ให้แน่ใจว่า n_tube2 เป็นจำนวนที่ลงตัวกับจำนวนคอลัมน์
        n_row2 = n_tube2 / column;

        % แสดงผลค่ากลางที่คำนวณได้
        fprintf('n_tube2 = %d, n_row2 = %d\n', n_tube2, n_row2);

        % คำนวณข้อมูลทางเศรษฐศาสตร์
        N_day = 360; % จำนวนวันที่ทำงานต่อปี (วัน/ปี)
        N = 20; % จำนวนปีที่ทำงาน (ปี)
        H_work = 10; % จำนวนชั่วโมงการทำงานต่อวัน (ชั่วโมง/วัน)
        i = 0.04; % อัตราราคาพลังงาน (0.xx)
        d = 0.016; % อัตราคิดลดของตลาด (0.xx)
        Ms = 0.06; % สัดส่วนของค่าใช้จ่ายประจำปีที่เกิดขึ้นในการดำเนินงานและบำรุงรักษาต่อค่าใช้จ่ายเดิม
        Rv = 0.1; % สัดส่วนของมูลค่าขายคืนต่อค่าใช้จ่ายเดิม
        cost_oil = 22; % ราคาน้ำมัน (บาท/ลิตร)
        CA = 3000; % ค่าใช้จ่ายแรกเริ่มที่ขึ้นอยู่กับพื้นที่ของ HPHE (บาท/ตารางเมตร)

        % คำนวณค่าพื้นฐาน
        H = N_day * H_work; % ชั่วโมง/ปี
        CE = cost_oil / 38581.5 * 3600 / 1000; % (บาท/Wh)
        A_HX = (pi * do * (Le + Lc)) * n_tube2; % พื้นที่ผิวของเครื่องแลกเปลี่ยนความร้อน

        % คำนวณค่า P1 และ P2
        P1 = compute_P1(N, i, d);
        P2 = 1 + P1 * Ms - Rv * (1 + d)^(-N);

        % แสดงผลการคำนวณค่าทางการเงิน
        fprintf('P1 = %.4f, P2 = %.4f, CE = %.4f\n', P1, P2, CE);

        % คำนวณค่าความต้านทาน
        [Z1, Z2, Z8, Z9, kx] = calculate_resistances(Le, Lc, do, di, n_tube, materials, heating_medium, cooling_medium, T_hi, T_ci, Q_h, Q_co, A_HX, column, Df_e, Df_c, characteristic_Of_HPHE, n_row2, ST, SL);

        Z_total = Z1 + Z9 + Z2 + Z8;

        % แสดงผลการคำนวณค่าความต้านทาน
        fprintf('Z1 = %.6f, Z2 = %.6f, Z8 = %.6f, Z9 = %.6f, Z_total = %.6f\n', Z1, Z2, Z8, Z9, Z_total);

        % คำนวณการถ่ายเทความร้อน
        [Q1, U, m_cp_min, Ch, Cc, C_st] = calculate_heat_transfer(Q_h, Q_co, Z_total, A_HX, T_hi, T_ci, heating_medium, cooling_medium);

        % แสดงผลการคำนวณค่าการถ่ายเทความร้อน
        fprintf('Q1 = %.2f, U = %.4f, m_cp_min = %.4f\n', Q1, U, m_cp_min);

        % คำนวณประสิทธิภาพ
        [Eff_c_opti, Eff_thermosy] = optimize_effectiveness(P1, P2, CE, H, m_cp_min, T_hi, T_ci, Cc, C_st, U, Q1, CA);

        % แสดงผลค่าต่างๆ ที่ได้คำนวณ
        fprintf('-------------------------------------\n');
        fprintf('การจัดเรียงท่อถูกปรับให้เข้ากับ %d tubes and %d rows.\n', n_tube2, n_row2);
        fprintf('พื้นที่การแลกเปลี่ยนความร้อน: %.2f m^2\n', A_HX);
        fprintf('ค่าการถ่ายเทความร้อนทั้งหมด: %.2f W/m^2K\n', U);
        fprintf('ประสิทธิภาพที่ได้รับการปรับปรุง: %.2f\n', Eff_c_opti);
        fprintf('ประสิทธิภาพของ thermosyphon: %.2f\n', Eff_thermosy);
        
    catch ME
        % การจัดการและแสดงข้อผิดพลาดที่เกิดขึ้น
        fprintf('เกิดข้อผิดพลาด: %s\n', ME.message);
        if strcmp(ME.identifier, 'MATLAB:dimagree')
            error('ขนาดของเมทริกซ์ต้องเท่ากัน: %s', ME.message);
        else
            rethrow(ME);
        end
    end
end


function P1 = compute_P1(N, i, d)
    try
        if i == d
            P1 = N / (1 + i);
        else
            P1 = (1 / (d - i)) * (1 - ((1 + i) / (1 + d))^N);
        end
        % Debug: Check P1 calculation
        fprintf('P1 calculation successful: P1 = %.4f\n', P1);
    catch
        error('Error in compute_P1 function');
    end
end

function [Z1, Z2, Z8, Z9, kx] = calculate_resistances(Le, Lc, do, di, n_tube, materials, heating_medium, cooling_medium, T_hi, T_ci, Q_h, Q_co, A_HX, column, Df_e, Df_c, characteristic_Of_HPHE, n_row2, ST, SL)
    try
        % Material properties
        kx = materials_tube(materials);
        fprintf('Material conductivity (kx) calculated: kx = %.4f\n', kx);
        
        % Calculate the resistances for evaporator and condenser sections
        fprintf('Calculating evaporator section resistance...\n');
        [Seo, heo] = calculate_section_resistance(Le, do, column, T_hi, Q_h, heating_medium, characteristic_Of_HPHE, n_row2, Df_e, ST, SL);
        fprintf('Evaporator section: Seo = %.4f, heo = %.4f\n', Seo, heo);
        
        fprintf('Calculating condenser section resistance...\n');
        [Sco, hco] = calculate_section_resistance(Lc, do, column, T_ci, Q_co, cooling_medium, characteristic_Of_HPHE, n_row2, Df_c, ST, SL);
        fprintf('Condenser section: Sco = %.4f, hco = %.4f\n', Sco, hco);
        
        Z1 = 1 / (heo * Seo);
        Z9 = 1 / (hco * Sco);
        Z2 = log(do / di) / (2 * n_tube * pi * Le * kx);
        Z8 = log(do / di) / (2 * n_tube * pi * Lc * kx);

        % Debug: Check individual resistances
        fprintf('Resistances calculated: Z1 = %.6f, Z2 = %.6f, Z8 = %.6f, Z9 = %.6f\n', Z1, Z2, Z8, Z9);

    catch ME
        fprintf('Error in calculate_resistances function at line %d: %s\n', ME.stack(1).line, ME.message);
        rethrow(ME);
    end
end


function [Seo, h] = calculate_section_resistance(L, do, column, T, Q, medium, characteristic_Of_HPHE, n_row2, Df, ST, SL)
    try
        % Debug: Start of section resistance calculation
        fprintf('Starting calculate_section_resistance...\n');
        fprintf('Inputs - L: %.4f, do: %.5f, column: %d, T: %.2f, Q: %.2f, medium: %d, characteristic_Of_HPHE: %d, n_row2: %d, Df: %.4f, ST: %.4f, SL: %.4f\n', ...
                L, do, column, T, Q, medium, characteristic_Of_HPHE, n_row2, Df, ST, SL);
        
        % Calculate cross-sectional flow area (A_cos_flow)
        if Df ~= 0
            A_cos_flow = (L + 0.01) * (do * column + (column * (Df - do)) + ((column - 1) * (0.47 - do - (Df - do)) + 0.05)); % m^2
        else
            A_cos_flow = (L + 0.01) * (do * column + ((column - 1) * (0.47 - do) + 0.05)); % m^2
        end
        
        % Debug: Check calculated flow area
        fprintf('Calculated A_cos_flow: %.6f m^2\n', A_cos_flow);

        % Calculate flow velocity (v)
        v = Q / (60 * 60 * A_cos_flow); % m/s
        fprintf('Calculated flow velocity (v): %.6f m/s\n', v);

        % Get properties of the medium
        if medium == 1
            [rhol, Cp, K, Pr, mev] = H2O_hp1(T);
        else
            [rhol, Cp, vr, K, alfa, Pr, Vo, mev] = AIR_hp(T);
        end
        
        % Debug: Check medium properties
        fprintf('Medium properties - rhol: %.4f, Cp: %.4f, K: %.6f, Pr: %.4f, mev: %.6f\n', rhol, Cp, K, Pr, mev);

        % Calculate Reynolds number (Re_D)
        Re_D = rhol * v * do / mev;
        fprintf('Calculated Reynolds number (Re_D): %.6f\n', Re_D);

        % Calculate Nusselt number (Nu_D) based on pipe alignment
        if characteristic_Of_HPHE == 1
            [c2] = Number_row_stagg(n_row2);
            [c1, m] = Reynold_num_Stagg(Re_D, ST, SL);
        else
            [c2] = Number_row_Alig(n_row2);
            [c1, m] = Reynold_num_Alig(Re_D);
        end
        
        % Debug: Check Nusselt number constants
        fprintf('Nusselt number constants - c1: %.4f, c2: %.4f, m: %.4f\n', c1, c2, m);

        % Final Nusselt number (Nu_D)
        Nu_D = c1 * c2 * Re_D^m * Pr^0.36 * (Pr / Pr)^0.25;
        fprintf('Calculated Nusselt number (Nu_D): %.6f\n', Nu_D);

        % Calculate effective surface area (Seo) and heat transfer coefficient (h)
        Seo = (pi * do * L) * column;
        h = Nu_D * K / do;
        
        % Debug: Final calculated values
        fprintf('Calculated Seo: %.6f m^2, h: %.6f W/m^2K\n', Seo, h);
        
    catch ME
        fprintf('Error in calculate_section_resistance function at line %d: %s\n', ME.stack(1).line, ME.message);
        rethrow(ME);
    end
end




function [Q1, U, m_cp_min, Ch, Cc, C_st] = calculate_heat_transfer(Q_h, Q_co, Z_total, A_HX, T_hi, T_ci, heating_medium, cooling_medium)
    try
        fprintf('Calculating heat transfer...\n');
        % Calculate heat transfer quantities
        dT = T_hi - T_ci;
        Q1 = dT / Z_total;
        U = 1 / (Z_total * A_HX);

        if heating_medium == 1
            [~, Cp_e] = H2O_hp1(T_hi);
        else
            [~, Cp_e] = AIR_hp(T_hi);
        end

        if cooling_medium == 1
            [~, Cp_c] = H2O_hp1(T_ci);
        else
            [~, Cp_c] = AIR_hp(T_ci);
        end

        m_ev = Q_h / 3600; % kg/s
        m_co = Q_co / 3600; % kg/s

        m_cp_min = min(m_ev * Cp_e, m_co * Cp_c);
        Ch = m_ev * Cp_e;
        Cc = m_co * Cp_c;
        C_st = Cc / Ch;

        % Debug: Check heat transfer calculations
        fprintf('Heat transfer calculated: Q1 = %.4f, U = %.4f, m_cp_min = %.4f\n', Q1, U, m_cp_min);
    catch
        error('Error in calculate_heat_transfer function');
    end
end

function [Eff_c_opti, Eff_thermosy] = optimize_effectiveness(P1, P2, CE, H, m_cp_min, T_hi, T_ci, Cc, C_st, U, Q1, CA)
    try
        fprintf('Starting optimize_effectiveness...\n');
        
        % Debug: Print inputs
        fprintf('Inputs - P1: %.4f, P2: %.4f, CE: %.4f, H: %.2f, m_cp_min: %.4f, T_hi: %.2f, T_ci: %.2f, Cc: %.4f, C_st: %.4f, U: %.4f, Q1: %.4f\n', ...
                P1, P2, CE, H, m_cp_min, T_hi, T_ci, Cc, C_st, U, Q1);
        
        % Calculate constants A and B
        A = P1 * CE * H * m_cp_min * (T_hi - T_ci);
        B = P2 * CA * Cc / ((1 - C_st) * U);

        % Debug: Print intermediate values
        fprintf('Intermediate values - A: %.4f, B: %.4f\n', A, B);

        % Check for potential division by zero or invalid operations
        if A == 0
            error('Division by zero: A is zero.');
        end
        if (1 - C_st) == 0
            error('Division by zero: (1 - C_st) is zero.');
        end

        % Optimized effectiveness calculation
        sqrt_term = sqrt((1 + C_st)^2 - 4 * C_st * (1 - B / A * (1 - C_st)));
        if sqrt_term < 0
            error('Invalid operation: square root of a negative number.');
        end
        
        Eff_c_opti = ((1 + C_st) - sqrt_term) / (2 * C_st);

        % Debug: Print optimized effectiveness
        fprintf('Optimized effectiveness (Eff_c_opti): %.4f\n', Eff_c_opti);

        % Calculate maximum possible heat transfer
        Q_max = m_cp_min * (T_hi - T_ci);
        if Q_max == 0
            error('Division by zero: Q_max is zero.');
        end
        
        % Calculate thermosyphon effectiveness
        Eff_thermosy = Q1 / Q_max;

        % Debug: Print thermosyphon effectiveness
        fprintf('Thermosyphon effectiveness (Eff_thermosy): %.4f\n', Eff_thermosy);

    catch ME
        fprintf('Error in optimize_effectiveness function at line %d: %s\n', ME.stack(1).line, ME.message);
        rethrow(ME);
    end
end


function [c2] = Number_row_stagg(n)
    % Check if n is less than 1
    if n < 1
        fprintf('Input is less than 1. Setting n to 1.\n');
        n = 1;
    % Check if n is greater than 9
    elseif n > 9
        fprintf('Input is greater than 9. Setting n to 9.\n');
        n = 9;
    end
    
    % c2 values for n from 1 to 9
    c2_values = [0.68, 0.75, 0.83, 0.89, 0.92, 0.95, 0.97, 0.98, 0.99];
    
    % Get the c2 value based on the (possibly adjusted) n
    c2 = c2_values(n);
    
    % Display the final c2 value
    fprintf('Using c2 value for n = %d: c2 = %.2f\n', n, c2);
end



function [c2] = Number_row_Alig(n)
    if n < 1
        error('No of tube')
    elseif n > 9
        error('Please input data (new column)')
    end
    c2 = [0.64 0.8 0.87 0.9 0.92 0.94 0.96 0.98 0.99];
    c2 = c2(n);
end


function [c1, m] = Reynold_num_Stagg(Re_D, ST, SL)
    Re = Re_D;
    if Re > 0
        if Re <= 100
            c1 = 0.90;
            m = 0.40;
        elseif Re <= 1000
            c1 = 0.51;
            m = 0.5;
        elseif Re <= 200000
            if (ST / SL) < 2
                c1 = 0.35 * (ST / SL)^(1/5);
            else
                c1 = 0.40;
            end
            m = 0.60;
        elseif Re <= 2000000
            c1 = 0.022;
            m = 0.84;
        end
    end
end


function [c1, m] = Reynold_num_Alig(Re_D)
    Re = Re_D;
    if Re > 10
        if Re <= 100
            c1 = 0.80;
            m = 0.40;
        elseif Re <= 1000
            c1 = 0.51;
            m = 0.5;
        elseif Re <= 200000
            c1 = 0.27;
            m = 0.63;
        elseif Re <= 2000000
            c1 = 0.021;
            m = 0.84;
        end
    end
end


function [kx] = materials_tube(materials)
% การคำนวณค่าการนำความร้อนของวัสดุที่ใช้ทำท่อ
if materials == 1
    kx = 16.2; % Steel
elseif materials == 2
    kx = 401; % Copper
else
    kx = 237; % Aluminium
end
end

function [n_row] = n_tube_column(n_row)
    if ~isnumeric(n_row) || isempty(n_row)
        error('Input must be a numeric value.');
    end

     disp(['Input number row value : ', num2str(n_row)]);
     
     if n_row < 1
        n_row = 1;
        disp('Input value is less than 1, setting to 1.');
     elseif n_row >= 1 && n_row <= 9
        n_row = ceil(n_row);
        disp('Input value exceeded 9, setting to 9.');
    elseif n_row > 9
        n_row = 9;
        disp('Input value exceeded 9, setting to 9.');
     
    end
end

function [rhol, Cp, vr, K, alfa, Pr, Vo, mev] = AIR_hp(T)
% ค่าอุณหภูมิและสมบัติต่างๆ ของอากาศที่อุณหภูมิ T
Temp_T = [100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750];
rhol_T = [3.5562, 2.3364, 1.7458, 1.3947, 1.1614, 0.9950, 0.8711, 0.7740, 0.6964, 0.6329, 0.5804, 0.5356, 0.4975, 0.4354];
Cp_T = [1.032, 1.012, 1.007, 1.006, 1.007, 1.009, 1.014, 1.021, 1.030, 1.040, 1.051, 1.063, 1.075, 1.087] * 1E3;
v_T = [2.00, 4.426, 7.590, 11.44, 15.89, 20.92, 26.41, 32.39, 38.79, 45.57, 52.69, 60.21, 68.1, 76.37] * 1E-6;
K_T = [9.34, 13.8, 18.1, 22.3, 26.3, 30.0, 33.8, 37.3, 40.7, 43.9, 46.9, 49.7, 52.4, 54.9] * 1E-3;
alfa_T = [2.54, 5.84, 10.3, 15.9, 22.5, 29.9, 38.3, 47.2, 56.7, 66.7, 76.9, 87.3, 98, 109] * 1E-6;
Pr_T = [0.786, 0.758, 0.737, 0.720, 0.707, 0.700, 0.690, 0.686, 0.684, 0.683, 0.685, 0.690, 0.695, 0.702];
mev_T = [71.1, 103.4, 132.5, 159.6, 184.6, 208.2, 230.1, 250.7, 270.1, 288.4, 305.8, 322.5, 338.8, 354.6] * 1E-7;

if T <= Temp_T(1)
    error('Working fluid is solid');
else
    i = 1;
    if T > Temp_T(end)
        i = length(Temp_T);
    else
        while T > Temp_T(i)
            i = i + 1;
        end
    end
end

T2 = Temp_T(i);
T1 = Temp_T(i - 1);
rhol = rhol_T(i - 1) + (rhol_T(i) - rhol_T(i - 1)) * (T - T1) / (T2 - T1);
Cp = Cp_T(i - 1) + (Cp_T(i) - Cp_T(i - 1)) * (T - T1) / (T2 - T1);
vr = v_T(i - 1) + (v_T(i) - v_T(i - 1)) * (T - T1) / (T2 - T1);
K = K_T(i - 1) + (K_T(i) - K_T(i - 1)) * (T - T1) / (T2 - T1);
alfa = alfa_T(i - 1) + (alfa_T(i) - alfa_T(i - 1)) * (T - T1) / (T2 - T1);
Pr = Pr_T(i - 1) + (Pr_T(i) - Pr_T(i - 1)) * (T - T1) / (T2 - T1);
mev = mev_T(i - 1) + (mev_T(i) - mev_T(i - 1)) * (T - T1) / (T2 - T1);
Vo = 1 / rhol;
end

function [rhol_l, Cp_l, k_l, mev_l, Pr_l] = H2O_hp1(T)

    Temp_T = [273.15 275 280 285 290 295 300 305 310 315 320 325 330 335 340 345 350 355 360 365 370 373.15 375 380 385 390 400 410 420 430 440 450 460 470 480 490];
    rhol_l_T = [1.000 1.000 1.000 1.000 1.000 1.000 0.999 0.999 0.998 0.998 0.998 0.997 0.997 0.996 0.996 0.996 0.995 0.995 0.994 0.994 0.993 0.992 0.991 0.990 0.989 0.988 0.983 0.978 0.973 0.968 0.963 0.958 0.953 0.948 0.943];
    Cp_l_T = [4.217 4.211 4.198 4.189 4.184 4.181 4.179 4.178 4.178 4.179 4.180 4.182 4.184 4.186 4.188 4.191 4.195 4.199 4.203 4.209 4.214 4.217 4.220 4.226 4.232 4.239 4.256 4.278 4.302 4.331 4.36 4.40 4.44 4.48 4.53] * 1E3;
    k_l_T = [569 574 582 590 598 606 613 620 628 634 640 645 650 656 660 668 668 671 674 677 679 680 681 683 685 686 688 688 688 685 682 678 673 667 660 651] * 1E-3;
    mev_l_T = [1750 1652 1422 1225 1080 959 855 769 695 631 577 528 489 453 420 389 365 343 324 306 289 279 274 260 248 237 217 200 185 173 162 152 143 136 129 124] * 1E-6;
    Pr_l_T = [13.5 12.9 12.2 11.6 11.0 10.4 9.9 9.4 8.9 8.4 8.0 7.6 7.2 6.8 6.5 6.2 5.9 5.6 5.3 5.1 4.8 4.6 4.4 4.2 4.0 3.8 3.4 3.1 2.8 2.6 2.4 2.2 2.0 1.8 1.7];
    
    if T <= Temp_T(1)
        error('Working fluid is solid')
    else
        i = 1;
        if T > Temp_T(end)
            i = length(Temp_T);
        else
            while T > Temp_T(i)
                i = i + 1;
            end
        end
    end

    T2 = Temp_T(i);
    T1 = Temp_T(i-1);
    rhol_l = rhol_l_T(i-1) + (rhol_l_T(i) - rhol_l_T(i-1)) * (T - T1) / (T2 - T1);
    Cp_l = Cp_l_T(i-1) + (Cp_l_T(i) - Cp_l_T(i-1)) * (T - T1) / (T2 - T1);
    k_l = k_l_T(i-1) + (k_l_T(i) - k_l_T(i-1)) * (T - T1) / (T2 - T1);
    mev_l = mev_l_T(i-1) + (mev_l_T(i) - mev_l_T(i-1)) * (T - T1) / (T2 - T1);
    Pr_l = Pr_l_T(i-1) + (Pr_l_T(i) - Pr_l_T(i-1)) * (T - T1) / (T2 - T1);
end


function [Ts] = H2O_hp_pressure(P)
% การคำนวณอุณหภูมิการอิ่มตัวของน้ำที่ความดัน P
Pres_T = [0.00611, 0.00697, 0.00990, 0.01387, 0.01917, 0.02617, 0.03531, 0.04712, 0.06221, 0.08132, 0.1053, 0.1351, 0.1719, 0.2167, 0.2713, 0.3372, 0.4163, 0.5100, 0.6209, 0.7514, 0.9040, 1.0133, 1.0815, 1.2869, 1.5233, 1.794, 2.455, 3.302, 4.370, 5.699, 7.333, 9.319, 11.71, 14.55, 17.90, 21.83];
Temp_T = [273.15, 275, 280, 285, 290, 295, 300, 305, 310, 315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 365, 370, 373.15, 375, 380, 385, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490];

if P <= Pres_T(1)
    error('Working fluid is solid');
else
    i = 1;
    if P > Pres_T(end)
        i = length(Pres_T);
    else
        while P > Pres_T(i)
            i = i + 1;
        end
    end
end

P2 = Pres_T(i);
P1 = Pres_T(i - 1);
Ts = (Temp_T(i - 1) + (Temp_T(i) - Temp_T(i - 1)) * (P - P1) / (P2 - P1)) - 273.15;
end

function [Pv, rhol_v, rhol_l, L, mev_l, sixmar_l, k_l, Cp_l] = R134a_hp(T)
% ค่าอุณหภูมิและสมบัติต่างๆ ของ R134a ที่อุณหภูมิ T
Temp_T = [273.15, 275.15, 277.15, 279.15, 281.15, 283.15, 285.15, 287.15, 289.15, 291.15, 293.15, 295.15, 297.15, 299.15, 301.15, 303.15, 305.15, 307.15, 309.15, 311.15, 313.15, 315.15, 317.15, 319.15, 321.15, 323.15, 325.15, 327.15, 329.15, 331.15, 333.15, 335.15, 337.15, 339.15, 341.15, 343.15, 345.15, 347.15, 349.15, 351.15, 353.15];
Pres_T = [2.928, 3.146, 3.377, 3.62, 3.876, 4.146, 4.43, 4.729, 5.043, 5.372, 5.717, 6.079, 6.458, 6.854, 7.269, 7.702, 8.154, 8.626, 9.118, 9.632, 10.17, 10.72, 11.3, 11.9, 12.53, 13.18, 13.85, 14.55, 15.28, 16.04, 16.82, 17.63, 18.47, 19.34, 20.24, 21.17, 22.13, 23.13, 24.16, 25.23, 26.33];
rhol_v_T = [14.43, 15.46, 16.56, 17.72, 18.94, 20.23, 21.58, 23.01, 24.52, 26.11, 27.78, 29.54, 31.39, 33.34, 35.38, 37.54, 39.8, 42.18, 44.68, 47.32, 50.09, 53, 56.06, 59.29, 62.69, 66.27, 70.05, 74.03, 78.24, 82.68, 87.38, 92.36, 97.64, 103.2, 109.2, 115.6, 122.4, 129.7, 137.5, 145.9, 155.1];
rhol_L_T = [1295, 1288, 1281, 1275, 1268, 1261, 1254, 1247, 1240, 1233, 1225, 1218, 1210, 1203, 1195, 1187, 1180, 1172, 1163, 1155, 1147, 1138, 1129, 1121, 1112, 1102, 1093, 1083, 1073, 1063, 1053, 1042, 1031, 1020, 1008, 996.2, 983.8, 970.8, 957.3, 943.1, 928.2];
h_L_T = [34.19, 36.88, 39.59, 42.3, 45.03, 47.77, 50.52, 53.28, 56.06, 58.85, 61.66, 64.48, 67.31, 70.16, 73.03, 75.91, 78.81, 81.73, 84.67, 87.62, 90.6, 93.6, 96.61, 99.66, 102.7, 105.8, 108.9, 112.1, 115.2, 118.5, 121.7, 125, 128.3, 131.6, 135, 138.5, 142, 145.5, 149.1, 152.8, 156.6];
h_V_T = [232.8, 234, 235.1, 236.3, 237.4, 238.5, 239.6, 240.7, 241.8, 242.9, 243.9, 245, 246, 247, 248, 249, 250, 250.9, 251.8, 252.7, 253.6, 254.5, 255.3, 256.1, 256.9, 257.6, 258.3, 259, 259.7, 260.3, 260.8, 261.3, 261.8, 262.2, 262.6, 262.8, 263.1, 263.2, 263.2, 263.2, 263];
mev_l_T = [287.4, 280.4, 273.6, 267, 260.6, 254.3, 248.3, 242.5, 236.8, 231.2, 225.8, 220.5, 215.4, 210.4, 205.5, 200.7, 196, 191.4, 186.9, 182.5, 178.2, 174, 169.8, 165.7, 161.7, 157.7, 153.8, 149.9, 146.1, 142.3, 138.6, 134.9, 131.2, 127.5, 123.9, 120.3, 116.7, 113.1, 109.4, 105.8, 102.1] * 1E-6;
Sur_Ten_T = [0.01156, 0.01127, 0.01099, 0.0107, 0.01042, 0.01014, 0.00986, 0.00958, 0.0093, 0.00903, 0.00876, 0.00848, 0.00821, 0.00795, 0.00768, 0.00742, 0.00715, 0.00689, 0.00664, 0.00638, 0.00613, 0.00588, 0.00563, 0.00538, 0.00513, 0.00489, 0.00465, 0.00441, 0.00418, 0.00395, 0.00372, 0.00349, 0.00327, 0.00305, 0.00283, 0.00261, 0.0024, 0.0022, 0.00199, 0.0018, 0.0016];
k_l_T = [0.09201, 0.09112, 0.09024, 0.08936, 0.08849, 0.08761, 0.08674, 0.08587, 0.08501, 0.08414, 0.08328, 0.08242, 0.08156, 0.0807, 0.07984, 0.07899, 0.07813, 0.07727, 0.07642, 0.07556, 0.07471, 0.07385, 0.073, 0.07214, 0.07128, 0.07042, 0.06956, 0.06869, 0.06783, 0.06696, 0.06609, 0.06521, 0.06433, 0.06345, 0.06256, 0.06167, 0.06078, 0.05988, 0.05898, 0.05807, 0.05717];
Cp_l_T = [1.341, 1.347, 1.352, 1.358, 1.364, 1.37, 1.377, 1.383, 1.39, 1.397, 1.405, 1.413, 1.421, 1.429, 1.437, 1.446, 1.456, 1.466, 1.476, 1.487, 1.498, 1.51, 1.523, 1.537, 1.551, 1.566, 1.582, 1.6, 1.618, 1.638, 1.66, 1.684, 1.71, 1.738, 1.769, 1.804, 1.843, 1.887, 1.938, 1.996, 2.065] * 1E3;

if T <= Temp_T(1)
    error('Working fluid is solid');
else
    i = 1;
    if T > Temp_T(end)
        i = length(Temp_T);
    else
        while T > Temp_T(i)
            i = i + 1;
        end
    end
end

T2 = Temp_T(i);
T1 = Temp_T(i - 1);
Pv = Pres_T(i - 1) + (Pres_T(i) - Pres_T(i - 1)) * (T - T1) / (T2 - T1);
rhol_v = rhol_v_T(i - 1) + (rhol_v_T(i) - rhol_v_T(i - 1)) * (T - T1) / (T2 - T1);
rhol_l = rhol_L_T(i - 1) + (rhol_L_T(i) - rhol_L_T(i - 1)) * (T - T1) / (T2 - T1);
L = h_V_T(i - 1) + (h_V_T(i) - h_V_T(i - 1)) * (T - T1) / (T2 - T1) - (h_L_T(i - 1) + (h_L_T(i) - h_L_T(i - 1)) * (T - T1) / (T2 - T1));
mev_l = mev_l_T(i - 1) + (mev_l_T(i) - mev_l_T(i - 1)) * (T - T1) / (T2 - T1);
sixmar_l = Sur_Ten_T(i - 1) + (Sur_Ten_T(i) - Sur_Ten_T(i - 1)) * (T - T1) / (T2 - T1);
k_l = k_l_T(i - 1) + (k_l_T(i) - k_l_T(i - 1)) * (T - T1) / (T2 - T1);
Cp_l = Cp_l_T(i - 1) + (Cp_l_T(i) - Cp_l_T(i - 1)) * (T - T1) / (T2 - T1);
end

function [Ts] = R134a_hp_pressure(P)
% การคำนวณอุณหภูมิการอิ่มตัวของ R134a ที่ความดัน P
Pres_T = [2.928, 3.146, 3.377, 3.62, 3.876, 4.146, 4.43, 4.729, 5.043, 5.372, 5.717, 6.079, 6.458, 6.854, 7.269, 7.702, 8.154, 8.626, 9.118, 9.632, 10.17, 10.72, 11.3, 11.9, 12.53, 13.18, 13.85, 14.55, 15.28, 16.04, 16.82, 17.63, 18.47, 19.34, 20.24, 21.17, 22.13, 23.13, 24.16, 25.23, 26.33];
Temp_T = [273.15, 275.15, 277.15, 279.15, 281.15, 283.15, 285.15, 287.15, 289.15, 291.15, 293.15, 295.15, 297.15, 299.15, 301.15, 303.15, 305.15, 307.15, 309.15, 311.15, 313.15, 315.15, 317.15, 319.15, 321.15, 323.15, 325.15, 327.15, 329.15, 331.15, 333.15, 335.15, 337.15, 339.15, 341.15, 343.15, 345.15, 347.15, 349.15, 351.15, 353.15];

if P <= Pres_T(1)
    error('Working fluid is solid');
else
    i = 1;
    if P > Pres_T(end)
        i = length(Pres_T);
    else
        while P > Pres_T(i)
            i = i + 1;
        end
    end
end

P2 = Pres_T(i);
P1 = Pres_T(i - 1);
Ts = Temp_T(i - 1) + (Temp_T(i) - Temp_T(i - 1)) * (P - P1) / (P2 - P1);
end



function [Pv, rhol_v, rhol_l, L, Cp_l, k_l, mev_l, sixmar_l] = H2O_hp(T)
    Temp_T = [273.15 275 280 285 290 295 300 305 310 315 320 325 330 335 340 345 350 355 360 365 370 373.15 375 380 385 390 400 410 420 430 440 450 460 470 480 490];
    Pres_T = [0.00611 0.00697 0.00990 0.01387 0.01917 0.02617 0.03531 0.04712 0.06221 0.08132 0.1053 0.1351 0.1719 0.2167 0.2713 0.3372 0.4163 0.5100 0.6209 0.7514 0.9040 1.0133 1.0815 1.2869 1.5233 1.794 2.455 3.302 4.370 5.699 7.333 9.319 11.71 14.55 17.90 21.83];
    rhol_v_T = [1/206.3 1/181.7 1/130.4 1/99.4 1/69.7 1/51.94 1/39.13 1/29.74 1/22.93 1/17.82 1/13.98 1/11.06 1/8.82 1/7.09 1/5.74 1/4.683 1/3.846 1/3.180 1/2.645 1/2.212 1/1.861 1/1.679 1/1.574 1/1.337 1/1.142 1/0.980 1/0.731 1/0.553 1/0.425 1/0.331 1/0.261 1/0.208 1/0.167 1/0.136 1/0.111 1/0.0922]; % kg/m^3
    rhol_l_T = [1 1 1 1 1 1 1/1.009 1/1.011 1/1.013 1/1.016 1/1.018 1/1.021 1/1.024 1/1.027 1/1.030 1/1.034 1/1.038 1/1.041 1/1.044 1/1.045 1/1.049 1/1.053 1/1.058 1/1.067 1/1.077 1/1.088 1/1.099 1/1.110 1/1.123 1/1.137 1/1.152 1/1.167 1/1.184]; % kg/m^3
    hfg_T = [2502 2497 2485 2473 2461 2449 2438 2426 2414 2402 2390 2378 2366 2354 2342 2329 2317 2304 2291 2278 2265 2257 2252 2239 2225 2212 2183 2153 2123 2091 2059 2024 1989 1951 1912 1870]; % kJ/kg
    mev_l_T = [1750 1652 1422 1225 1080 959 855 769 695 631 577 528 489 453 420 389 365 343 324 306 289 279 274 260 248 237 217 200 185 173 162 152 143 136 129 124]; % mPa·s
    sixmar_l_T = [0.075 0.074 0.073 0.071 0.070 0.069 0.067 0.066 0.064 0.063 0.061 0.060 0.059 0.057 0.056 0.055 0.054 0.053 0.052 0.051 0.050 0.049 0.048 0.047 0.046 0.044 0.043 0.042 0.041 0.040 0.039 0.038 0.037 0.036]; % N/m
    k_l_T = [569 574 582 590 598 606 613 620 628 634 640 645 650 656 660 668 668 671 674 677 679 680 681 683 685 686 688 688 688 685 682 678 673 667 660 651]; % W/m·K
    Cp_l_T = [4.217 4.211 4.198 4.189 4.184 4.181 4.179 4.178 4.178 4.179 4.180 4.182 4.184 4.186 4.188 4.191 4.195 4.199 4.203 4.209 4.214 4.217 4.220 4.226 4.232 4.239 4.256 4.278 4.302 4.331 4.36 4.40 4.44 4.48 4.53]; % kJ/kg·K

     % แสดงค่าที่ป้อนเข้า
    fprintf('Input Temperature T: %f\n', T);
    
    % ตรวจสอบค่า T และกำหนดค่า i
    if T <= Temp_T(1)
        error('Working fluid is solid');
    elseif T >= Temp_T(end)
        i = length(Temp_T);
    else
        i = find(Temp_T >= T, 1);
    end
    
    % แสดงค่าของ i
    fprintf('Initial index i: %d\n', i);
    
    % กำหนดค่า i ให้อยู่ในขอบเขตที่ถูกต้อง
    i = min(max(i, 2), length(Temp_T));
    i_minus1 = max(i - 1, 1);
    
    % แสดงค่าของ i และ i-1 หลังการปรับ
    fprintf('Adjusted index i: %d\n', i);
    fprintf('Adjusted index i-1: %d\n', i_minus1);
    
    % ตรวจสอบขอบเขตของ i และ i_minus1 สำหรับอาเรย์ที่เกี่ยวข้อง
    if i > length(Pres_T)
        error('Index i is out of bounds for Pres_T. i = %d, Length of Pres_T = %d', i, length(Pres_T));
    end
    if i_minus1 < 1 || i_minus1 > length(Pres_T)
        error('Index i-1 is out of bounds for Pres_T. i-1 = %d, Length of Pres_T = %d', i_minus1, length(Pres_T));
    end
    if i > length(rhol_v_T)
        error('Index i is out of bounds for rhol_v_T. i = %d, Length of rhol_v_T = %d', i, length(rhol_v_T));
    end
    if i_minus1 < 1 || i_minus1 > length(rhol_v_T)
        error('Index i-1 is out of bounds for rhol_v_T. i-1 = %d, Length of rhol_v_T = %d', i_minus1, length(rhol_v_T));
    end
    if i > length(rhol_l_T)
        error('Index i is out of bounds for rhol_l_T. i = %d, Length of rhol_l_T = %d', i, length(rhol_l_T));
    end
    if i_minus1 < 1 || i_minus1 > length(rhol_l_T)
        error('Index i-1 is out of bounds for rhol_l_T. i-1 = %d, Length of rhol_l_T = %d', i_minus1, length(rhol_l_T));
    end
    if i > length(hfg_T)
        error('Index i is out of bounds for hfg_T. i = %d, Length of hfg_T = %d', i, length(hfg_T));
    end
    if i_minus1 < 1 || i_minus1 > length(hfg_T)
        error('Index i-1 is out of bounds for hfg_T. i-1 = %d, Length of hfg_T = %d', i_minus1, length(hfg_T));
    end
    if i > length(mev_l_T)
        error('Index i is out of bounds for mev_l_T. i = %d, Length of mev_l_T = %d', i, length(mev_l_T));
    end
    if i_minus1 < 1 || i_minus1 > length(mev_l_T)
        error('Index i-1 is out of bounds for mev_l_T. i-1 = %d, Length of mev_l_T = %d', i_minus1, length(mev_l_T));
    end
    if i > length(sixmar_l_T)
        error('Index i is out of bounds for sixmar_l_T. i = %d, Length of sixmar_l_T = %d', i, length(sixmar_l_T));
    end
    if i_minus1 < 1 || i_minus1 > length(sixmar_l_T)
        error('Index i-1 is out of bounds for sixmar_l_T. i-1 = %d, Length of sixmar_l_T = %d', i_minus1, length(sixmar_l_T));
    end
    if i > length(k_l_T)
        error('Index i is out of bounds for k_l_T. i = %d, Length of k_l_T = %d', i, length(k_l_T));
    end
    if i_minus1 < 1 || i_minus1 > length(k_l_T)
        error('Index i-1 is out of bounds for k_l_T. i-1 = %d, Length of k_l_T = %d', i_minus1, length(k_l_T));
    end
    if i > length(Cp_l_T)
        error('Index i is out of bounds for Cp_l_T. i = %d, Length of Cp_l_T = %d', i, length(Cp_l_T));
    end
    if i_minus1 < 1 || i_minus1 > length(Cp_l_T)
        error('Index i-1 is out of bounds for Cp_l_T. i-1 = %d, Length of Cp_l_T = %d', i_minus1, length(Cp_l_T));
    end
    
    % แสดงการตรวจสอบขอบเขต
    fprintf('Index i is within bounds for all relevant arrays.\n');
    
    % คำนวณค่าต่าง ๆ
    T2 = Temp_T(i);
    T1 = Temp_T(i_minus1);

    % แสดงค่าที่ใช้ในการคำนวณ
    fprintf('Temperature Range: T1 = %f, T2 = %f\n', T1, T2);

    Pv = Pres_T(i-1) + (Pres_T(i) - Pres_T(i-1)) * (T - T1) / (T2 - T1);
    rhol_v = rhol_v_T(i-1) + (rhol_v_T(i) - rhol_v_T(i-1)) * (T - T1) / (T2 - T1);
    rhol_l = rhol_l_T(i-1) + (rhol_l_T(i) - rhol_l_T(i-1)) * (T - T1) / (T2 - T1);
    L = hfg_T(i-1) + (hfg_T(i) - hfg_T(i-1)) * (T - T1) / (T2 - T1);
    mev_l = mev_l_T(i-1) + (mev_l_T(i) - mev_l_T(i-1)) * (T - T1) / (T2 - T1);
    sixmar_l = sixmar_l_T(i-1) + (sixmar_l_T(i) - sixmar_l_T(i-1)) * (T - T1) / (T2 - T1);
    k_l = k_l_T(i-1) + (k_l_T(i) - k_l_T(i-1)) * (T - T1) / (T2 - T1);
    Cp_l = Cp_l_T(i-1) + (Cp_l_T(i) - Cp_l_T(i-1)) * (T - T1) / (T2 - T1);

    % แสดงค่าผลลัพธ์ที่คำนวณได้
    fprintf('Calculated Values:\n');
    fprintf('Pv: %f\n', Pv);
    fprintf('rhol_v: %f\n', rhol_v);
    fprintf('rhol_l: %f\n', rhol_l);
    fprintf('L: %f\n', L);
    fprintf('mev_l: %f\n', mev_l);
    fprintf('sixmar_l: %f\n', sixmar_l);
    fprintf('k_l: %f\n', k_l);
    fprintf('Cp_l: %f\n', Cp_l);
end