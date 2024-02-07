function [Polar] = Airfoil_Runner(NACA_Number,test_re,test_mach,points)
try
        xf = XFOIL;
        xf.KeepFiles = false; % Set it to true to keep all intermediate files created (Airfoil, Polars, ...)
        xf.Visible = false;    % Set it to false to hide XFOIL plotting window
        xf.Airfoil =  Airfoil.createNACA4(NACA_Number,150);
        %Add five filtering steps to smooth the airfoil coordinates and help convergence
        xf.addFiltering(5);
        %Switch to OPER mode, and set Reynolds and mach
        xf.addOperation(test_re, test_mach);
        %Set maximum number of iterations
        xf.addIter(300)
        %Initializate the calculations
        xf.addAlpha(0,true);
        %Create a new polar
        xf.addPolarFile('Polar.txt');
        % %Calculate a sequence of angle of attack, from 0 to 25 degrees, step size of 0.1 degrees
        xf.addAlpha(points);
        %Close the polar file
        xf.addClosePolarFile;
        %And finally add the action to quit XFOIL
        xf.addQuit;
        xf.run
        %fprintf('Running NACA%s\n',NACA_Number)
        finished = xf.wait(10);
    
        try
            fid = fopen(char(xf.PolarFiles));
            status = fseek(fid, 429, 'bof');
            data=textscan(fid,'%f%f%f%f%f%f%f');
            fclose(fid);
            if not(isempty(data)) && status == 0
                data=cell2mat(data);
            else
                data = [];
            end
    
            nn=any(isnan(data),2);
            data(nn,:) = [];
            Polar.Alpha   = data(:,1);
            Polar.CL      = data(:,2);
            Polar.CD      = data(:,3);
            xf.Polars = Polar;
        catch
  

        end

        [status, cmdout] = system('tasklist /FI "IMAGENAME eq xfoil.exe"');
        pattern = 'xfoil.exe\s+(\d+)';
        matches = regexp(cmdout, pattern, 'tokens');
        pid = str2double(matches{1}{1});

        % Now kill the process using its PID
        if ~isempty(pid)
            system(['taskkill /F /PID ' num2str(pid)]);
            fprintf('Killed Xfoil Process\n')
        end
        delete(xf.ActionsFile)
        delete(char(xf.AirfoilFile))
        delete(char(xf.PolarFiles))

        
catch
        delete(char(xf.PolarFiles))
        delete(char(xf.ActionsFile))
        delete(char(xf.AirfoilFile))

end

