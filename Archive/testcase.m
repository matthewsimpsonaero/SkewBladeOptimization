xf = XFOIL;
xf.KeepFiles = false; % Set it to true to keep all intermediate files created (Airfoil, Polars, ...)
xf.Visible = false;    % Set it to false to hide XFOIL plotting window
xf.Airfoil =  Airfoil.createNACA4('0006',150);
%Add five filtering steps to smooth the airfoil coordinates and help convergence
xf.addFiltering(5);

%Switch to OPER mode, and set Reynolds = 3E7, Mach = 0.1
xf.addOperation(10000, 0.1);
%Set maximum number of iterations
xf.addIter(100)
%Initializate the calculations
xf.addAlpha(0,true);
%Create a new polar
xf.addPolarFile('Polar.txt');
% %Calculate a sequence of angle of attack, from 0 to 25 degrees, step size of 0.1 degrees
xf.addAlpha(0:2.5:25);
%Close the polar file
xf.addClosePolarFile;
%And finally add the action to quit XFOIL
xf.addQuit;

%% Now we're ready to run XFOIL
xf.run
disp('Running XFOIL, please wait...')

%% Wait up to 100 seconds for it to finish... 
%It is possible to run more than one XFOIL instance at the same time
finished = xf.wait(100); 

%% If successfull, read and plot the polar
if finished
    disp('XFOIL analysis finished.')
end

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
            
            n=any(isnan(data),2);
            data(n,:) = [];
            Polar.Alpha   = data(:,1);
            Polar.CL      = data(:,2);
            Polar.CD      = data(:,3);
            Polar.CDp     = data(:,4);
            Polar.CM      = data(:,5);
            Polar.Top_Xtr = data(:,6);
            Polar.Bot_Xtr = data(:,7);
            
            xf.Polars = Polar;
        catch
  
end
delete(char(xf.PolarFiles))
delete(char(xf.ActionsFile))
delete(char(xf.AirfoilFile))