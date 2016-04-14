function ScanMatchInfo = ScanMatch_Struct(ScanMatchInfo)
% SCANMATCH_STRUCT is an Interactive GUI which prepares a structure to be 
% used with the ScanMatch toolbox. This function creates a substitution 
% matrix based on the Euclidian distance between each RoI and a RoI grid
% mask.  
% 
% This toolbox will output a structure with the following fields:
%
% ScanMatchInfo = 
% 
%           Xres: 1024
%           Yres: 768
%           Xbin: 12
%           Ybin: 8
%     RoiModulus: 12
%      Threshold: 3.5000
%       GapValue: 0
%        TempBin: 0
%      SubMatrix: [96x96 double]
%           mask: [768x1024 double]
%
% 
%   Part of the ScanMatch toolbox
%   Written by Filipe Cristino 
%   $Version: 1.00 $  $Date: 10/09/2009

% If ScanMatchInfo is inputed then display its parameters
if nargin == 1 
    % Check if the ScanMatchInfo structure is of the right form
    Ok = ScanMatch_CheckStructure(ScanMatchInfo);
else
    ScanMatchInfo.Xres = 1024;
    ScanMatchInfo.Yres = 768;
    ScanMatchInfo.Xbin = 12.0;
    ScanMatchInfo.Ybin = 8.0;
    ScanMatchInfo.RoiModulus = ScanMatchInfo.Xbin;
    ScanMatchInfo.Threshold = 3.5;
    ScanMatchInfo.GapValue = 0;
    ScanMatchInfo.TempBin = 0;
    % Compute substitution matrix
    ScanMatchInfo.SubMatrix = ScanMatch_CreateSubMatrix(ScanMatchInfo.Xbin,...
            ScanMatchInfo.Ybin, ScanMatchInfo.Threshold);        
end


% Create and setup figure
fh = figure('Position',[200 200 650 550],'Toolbar','none', ...
            'name', 'ScanMatch, Info structure tool','NumberTitle', 'off',...
            'MenuBar', 'none', 'Color', [0.95,0.95,0.95],'wvisual','27');


ah1 = axes('Parent',fh,'units','pixels',...
           'Position',[20 140 500 400 ],...
           'Visible', 'off',...
           'DataAspectRatioMode','auto'...
       );  
   

% Display mask
display_mask()



% Create the menus
% Create menu to swap display 
hpD = uibuttongroup('Title','Display:','FontSize',12,...
             'units','pixels',...
             'tag', 'display',...
             'Position',[525 340 120 80 ]);
         
DispRoi = uicontrol('Parent', hpD, 'Style', 'radiobutton' ,...
               'String', 'RoI mask',...
               'Tag', 'disp_RoiMask', ...
               'Position', [0 30 116 20]...
               );

DispSub = uicontrol('Parent', hpD, 'Style', 'radiobutton' ,...
               'String', 'Substitution Matrix',...
               'Tag', 'disp_SubMatrix', ...
               'Position', [0 5 116 20]...
               );           
set(hpD,'SelectionChangeFcn',@popup_menu_Callback);

% Create panel to hold the parameters
hp = uipanel('Title','Parameters','FontSize',12,...
             'units','pixels',...
             ...'BackgroundColor','white',...
             'Position',[0 0 650 120 ]);

% Create panel and controls for the experiment size        
hpR = uipanel('Parent', hp,'Title','Stimuli resolution','FontSize',10,...
             'units','pixels',...
             'BorderType', 'none', ...
             'Position',[40 55 190 40 ]);
         
uicontrol('Parent', hpR, 'Style','text',...
            'String','X = ',...
            'Position', [0 0 20 20]);

        
XresUI = uicontrol('Parent', hpR, 'Style', 'edit' ,...
               'String', ScanMatchInfo.Xres,...
               'Tag', 'Xres', ...
               'Position', [25 0 50 20],...
               'Callback',{@popup_menu_Callback} ...
               );

uicontrol('Parent', hpR, 'Style','text',...
            'String','Y = ',...
            'Position', [80 0 20 20]);

        
YresUI = uicontrol('Parent', hpR, 'Style', 'edit' ,...
               'String', ScanMatchInfo.Yres,...
               'Tag', 'Yres', ...
               'Position', [105 0 50 20],...
               'Callback',{@popup_menu_Callback} ...
               );
           
uicontrol('Parent', hpR, 'Style','text',...
            'String','pixels ',...
            'Position', [156 0 40 20]);
           
% Create panel and controls for the Number of bins 
hpR = uipanel('Parent', hp,'Title','Number of bins','FontSize',10,...
             'units','pixels',...
             ...'BackgroundColor','white',...
             'BorderType', 'none', ...
             'Position',[40 10 190 40 ]);
         
uicontrol('Parent', hpR, 'Style','text',...
            'String','X = ',...
            'Position', [0 0 20 20]);

        
XbinUI = uicontrol('Parent', hpR, 'Style', 'edit' ,...
               'String', ScanMatchInfo.Xbin,...
               'Tag', 'Xbin', ...
               'Position', [25 0 50 20],...
               'Callback',{@popup_menu_Callback} ...
               );

uicontrol('Parent', hpR, 'Style','text',...
            'String','Y = ',...
            'Position', [80 0 20 20]);

        
YbinUI = uicontrol('Parent', hpR, 'Style', 'edit' ,...
               'String', ScanMatchInfo.Ybin,...
               'Tag', 'Ybin', ...
               'Position', [105 0 50 20],...
               'Callback',{@popup_menu_Callback} ...
               );
           
% Create controls for the Substitution Matrix Threshold 
uicontrol('Parent', hp, 'Style','text',...
            'String','Substitution Matrix Threshold = ',...
            'HorizontalAlignment', 'right', ...
            'Position', [230 80 150 20]);

        
ThresholdUI = uicontrol('Parent', hp, 'Style', 'edit' ,...
               'String', ScanMatchInfo.Threshold,...
               'Tag', 'Threshold', ...
               'Position', [390 80 50 20],...
               'Callback',{@popup_menu_Callback} ...
               );
           
% Create controls for the Gap Value 
uicontrol('Parent', hp, 'Style','text',...
            'HorizontalAlignment', 'right', ...
            'String','Gap Value = ',...
            'Position', [475 80 100 20]);

        
GapValueUI = uicontrol('Parent', hp, 'Style', 'edit' ,...
               'String', ScanMatchInfo.GapValue,...
               'Tag', 'GapValue', ...
               'Position', [590 80 50 20],...
               'Callback',{@popup_menu_Callback} ...
               );
           
% Create controls for the Temporal binning 
uicontrol('Parent', hp, 'Style','text',...
            'String','Temporal bin size (ms) = ',...
            'HorizontalAlignment', 'right', ...
            'Position', [460 55 120 20]);

        
TempBinUI = uicontrol('Parent', hp, 'Style', 'edit' ,...
               'String', ScanMatchInfo.TempBin,...
               'Tag', 'TempBin', ...
               'Position', [590 55 50 20],...
               'Callback',{@popup_menu_Callback} ...
               );
           
% Create controls for the RoiModulus       
uicontrol('Parent', hp, 'Style','text',...
            'HorizontalAlignment', 'right', ...
            'String','RoI Modulus = ',...
            'Position', [230 55 150 20]);

        
RoiModulusUI = uicontrol('Parent', hp, 'Style', 'edit' ,...
               'String', ScanMatchInfo.RoiModulus,...
               'Tag', 'RoiModulus', ...
               'Position', [ 390 55 50 20],...
               'Callback',{@popup_menu_Callback} ...
               );
           
% Create button to exit 
uicontrol('Parent', hp, 'Style', 'pushbutton' ,...
               'String', 'Save and exit', ...
               'Tag', 'Exit', ...
               'Position', [ 350 15 200 30],...
               'Callback',{@popup_menu_Callback} ...
               );            
           
    function popup_menu_Callback(source,eventdata)
        % Nested function hence all the variables are globals.
        
        pop_menu = get(source,'Tag');
        switch pop_menu
            case ('Xres')
                Xres_text = get(source,'String');
                Xres_new = str2double(Xres_text);
                if (Xres_new > 1 && Xres_new < 3000)
                    ScanMatchInfo.Xres = Xres_new;
                    display_mask();
                else
                    msgbox('Wrong value parsed', 'error', 'error')
                    set(XresUI, 'string', ScanMatchInfo.Xres)
                end
            case ('Yres')
                Yres_text = get(source,'String');
                Yres_new = str2double(Yres_text);
                if (Yres_new > 1 && Yres_new < 3000)
                    ScanMatchInfo.Yres = Yres_new;
                    display_mask();
                else
                    msgbox('Wrong value parsed', 'error', 'error')
                    set(YresUI, 'string', ScanMatchInfo.Yres)
                end
            case ('Xbin')
                Xbin_text = get(source,'String');
                Xbin_new = str2double(Xbin_text);
                if (Xbin_new > 1 && Xbin_new < 3000 && ((Xbin_new *  ScanMatchInfo.Ybin) < 677)) 
                    ScanMatchInfo.Xbin = Xbin_new;
                    ScanMatchInfo.SubMatrix = ScanMatch_CreateSubMatrix(ScanMatchInfo.Xbin,...
                        ScanMatchInfo.Ybin, ScanMatchInfo.Threshold);
                    display_mask();
                else
                    if ((Xbin_new *  ScanMatchInfo.Ybin) > 676)
                        msgbox('Number of RoIs cannot be superior to 676!', 'error', 'error')
                    else
                        msgbox('Wrong value parsed', 'error', 'error')
                    end
                    set(XbinUI, 'string', ScanMatchInfo.Xbin)
                end
            case ('Ybin')
                Ybin_text = get(source,'String');
                Ybin_new = str2double(Ybin_text);
                if (Ybin_new > 1 && Ybin_new < 3000 && ((Ybin_new *  ScanMatchInfo.Xbin) < 677))
                    ScanMatchInfo.Ybin = Ybin_new;
                     ScanMatchInfo.SubMatrix = ScanMatch_CreateSubMatrix(ScanMatchInfo.Xbin,...
                        ScanMatchInfo.Ybin, ScanMatchInfo.Threshold);
                    display_mask();
                    
                else
                     if ((Ybin_new *  ScanMatchInfo.Xbin) > 676)
                        msgbox('Number of RoIs cannot be superior to 676!', 'error', 'error')
                    else
                        msgbox('Wrong value parsed', 'error', 'error')
                    end
                    set(YbinUI, 'string', ScanMatchInfo.Ybin)
                end
            case ('Threshold')
                Threshold_text = get(source,'String');
                Threshold_new = str2double(Threshold_text);
                if (Threshold_new > 1 && Threshold_new < 100)
                    ScanMatchInfo.Threshold = Threshold_new;
                    % Compute substitution matrix
                    ScanMatchInfo.SubMatrix = ScanMatch_CreateSubMatrix(ScanMatchInfo.Xbin,...
                        ScanMatchInfo.Ybin, ScanMatchInfo.Threshold);
                    % display the substitution matrix if onScreen
                    if get(DispSub, 'value') == 1
                        display_SubMatrix()
                                
                    end
                else
                    msgbox('Wrong value parsed', 'error', 'error')
                    set(ThresholdUI, 'string', ScanMatchInfo.Threshold)
                end
            case ('RoiModulus')
                RoiModulus_text = get(source,'String');
                RoiModulus_new = str2double(RoiModulus_text);
                if (RoiModulus_new > 1 && RoiModulus_new < 27)
                    ScanMatchInfo.RoiModulus = RoiModulus_new;
                    display_mask();
                else
                    msgbox('Wrong value parsed', 'error', 'error')
                    set(RoiModulusUI, 'string', ScanMatchInfo.RoiModulus)
                end
            case('GapValue')
                GapValue_text = get(source,'String');
                GapValue_new = str2double(GapValue_text);
                if (GapValue_new > -100 && GapValue_new < 100)
                    ScanMatchInfo.GapValue = GapValue_new;
                    
                else
                    msgbox('Wrong value parsed', 'error', 'error')
                    set(GapValueUI, 'string', ScanMatchInfo.GapValue)
                end
            case('TempBin')
                TempBin_text = get(source,'String');
                TempBin_new = str2double(TempBin_text);
                if (TempBin_new >= 0 && TempBin_new < 1000)
                    ScanMatchInfo.TempBin = TempBin_new;
                    
                else
                    msgbox('Wrong value parsed', 'error', 'error')
                    set(TempBinUI, 'string', ScanMatchInfo.TempBin)
                end    
            case ('display')
                switch get(eventdata.NewValue,'Tag')
                    case ('disp_SubMatrix')
                        display_SubMatrix()
                    case ('disp_RoiMask')
                        display_mask()
                end
            case ('Exit')
                disp('The following ScanMatchInfo structure has been created:')
                disp(ScanMatchInfo);
                assignin('base', 'ScanMatchInfo', ScanMatchInfo);
                close(fh)
        end
    end

    function display_mask()
        %Create mask output   
        ScanMatchInfo.mask = ScanMatch_GridMask(ScanMatchInfo.Xres,ScanMatchInfo.Yres, ...
            ScanMatchInfo.Xbin, ScanMatchInfo.Ybin);
        imagesc(ScanMatchInfo.mask);
        axis off
                       
        %Annotate Rois
        Y_ratio = ScanMatchInfo.Yres/ScanMatchInfo.Ybin;
        X_ratio = ScanMatchInfo.Xres/ScanMatchInfo.Xbin;
        ind = 1;
        for i= Y_ratio/2:Y_ratio:ScanMatchInfo.Yres
            for j= X_ratio/2:X_ratio:ScanMatchInfo.Xres
                st = ScanMatch_NumToDoubleStr(ind,ScanMatchInfo.RoiModulus);
                text(j,i, st);
                ind = ind+1;
            end
        end
    end   

    function display_SubMatrix()
        %Create mask output   
        
        imagesc(ScanMatchInfo.SubMatrix) 
        max(ScanMatchInfo.SubMatrix(:))
        colorbar
        axis off
        
        
    end        

end