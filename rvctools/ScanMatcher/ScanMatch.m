function [score alignment] = ScanMatch(seq1, seq2, ScanMatchInfo, varargin)
%SCANMATCH compares two sequences using the Needleman-Wunsch 
% algorithm providing a GUI is needed to display the results. 
%   SCORE = SCANMATCH(SEQ1, SEQ2, SCANMATCHINFO) where SEQ1 and SEQ2 are
%   two double strings (where a single fixation is represented by two
%   characters). SCANMATCHINFO is a structure containing all the parameters
%   needed to perform the comparison and must contain the following fields:
%
%           .Xres
%           .Yres
%           .Xbin
%           .Ybin
%     .RoiModulus
%      .Threshold
%       .GapValue
%        .TempBin
%      .SubMatrix
%           .mask
%
%   As an extra input parameter 'ShowViewer'm can be set to 1 to display the 
%   GUI. If not set at all or set to zero the Gui will be not shown.
%   
%   SCORE is the normalised result of the alignment. The closest SCORE is
%   to one the more related the two sequences are. Alignment output the
%   alignment of the two sequences. This is quite computationally expensive
%   hence if not needed call the function as follow:
%           score = ScanMatch(seq1, seq2, ScanMatchInfo)
%
%   Part of the ScanMatch toolbox
%   Written by Filipe Cristino 
%   $Version: 1.00 $  $Date: 10/09/2009



% Check the input parameters
okargs = {'ShowViewer'};
if nargin == 3
    ShowViewer = 0; % if not specified then don't show
else
    if strcmp(varargin(1), okargs)
        ShowViewer = varargin{2};
    else
        error('ScanMatch:InvalidInput','Invalid parameter: %s', char(varargin(1)));
    end
end

% Check if the ScanMatchInfo structure is of the right form
Ok = ScanMatch_CheckStructure(ScanMatchInfo);
   
% use numerical arrays for easy indexing
if ischar(seq1)
    intseq1 = uint8(ScanMatch_DoubleStrToNum(seq1,ScanMatchInfo.RoiModulus));
else
    intseq1 = double(seq1);
    seq1 = ScanMatch_NumToDoubleStr(intseq1,ScanMatchInfo.RoiModulus);
end
if ischar(seq2)
    intseq2 = uint8(ScanMatch_DoubleStrToNum(seq2,ScanMatchInfo.RoiModulus));
else
    intseq2 = double(seq2);
    seq2 = ScanMatch_NumToDoubleStr(intseq2,ScanMatchInfo.RoiModulus);
end

% Check for the sequences lenghts >0 and that they are even 

m = length(seq1);
n = length(seq2);

if  ~n||~m || mod(m,2)||mod(n,2)
    error('ScanMatch:InvalidLengthSequences','Length of input sequences must be greater than 0 and even');
end

m = length(intseq1);
n = length(intseq2);

% Perform the Needleman-Wunsch alogorithm
[score path step F] = ScanMatch_nwAlgo(intseq1,m,intseq2,n,ScanMatchInfo.SubMatrix,ScanMatchInfo.GapValue);

% Normalise output score
max_sub = max(ScanMatchInfo.SubMatrix(:));
scale = max_sub * max(m,n);
score = score / scale; 

% if alignment output is not needed then exit now (faster)
if nargout<=1 && ~ShowViewer
    return
end

path(step:end,:) = []; % step-1 is the length of the alignment
path = flipud(path);

%double the path to take into account the double string
p1 = interp1(1:step-1,double(path(:,1)),fix(1:0.5:step-1+0.5),'nearest');
p2 = interp1(1:step-1,double(path(:,2)),fix(1:0.5:step-1+0.5),'nearest');
path = [p1' p2'];

% setting the size of the alignment
n_l = (step-1)*2;
alignment = repmat(('- -')',1,n_l);

% adding sequence to alignment
alignment(1,path(:,1)>0) = seq1;
alignment(3,path(:,2)>0) = seq2;

% find locations where there are no gaps
h=find(all(path>0,2));
noGaps1=ScanMatch_DoubleStrToNum(alignment(1,h),ScanMatchInfo.RoiModulus);
noGaps2=ScanMatch_DoubleStrToNum(alignment(3,h),ScanMatchInfo.RoiModulus);

% score pairs with no gap
value = ScanMatchInfo.SubMatrix(sub2ind(size(ScanMatchInfo.SubMatrix),double(noGaps1),double(noGaps2)));
value = interp1(1:length(value),value,fix(1:0.5:length(value)+0.5),'nearest');

% insert symbols of the match string into the alignment
alignment(2,h(value>=0)) = ':';
alignment(2,h(value==max(ScanMatchInfo.SubMatrix(:)))) = '|';

if ShowViewer == 1 % display the Gui viewer
    
    % Set some variable
    res_size = [800 380];

    % Create and setup figure
    fh = figure('Position',[200 200 res_size],'Toolbar','none', ...
                'name', 'String Alignment Viewer','NumberTitle', 'off',...
                'MenuBar', 'none', 'Color', [0.95,0.95,0.95],'wvisual','27');

    uicontrol('Style','text', ...
                    'String','Alignment Results',...
                    'Position', [res_size(1)/2 - 100 res_size(2)-25 200 25], 'Fontsize', 16);

    ah2 = axes('Parent',fh,'units','pixels',...
               'Position',[40 40 280 180 ],...
               'Visible', 'off',...
               'DataAspectRatioMode','auto'...
           );               

    
    % Display the Score matrix and winning path   
    F=scale.*max(F(2:end,2:end,:),[],3);
    clim=max(max(max(abs(F(~isinf(F))))),eps);
    imagesc(F,[-clim clim]);
    colormap(privateColorMap(1));
    set(colorbar,'YLim',[min([F(:);-eps]) max([F(:);eps])])
    title('Scoring Space and Winning Path')
    xlabel('Sequence 1')
    ylabel('Sequence 2')
    hold on
    plot(path(all(path>0,2),1),path(all(path>0,2),2),'k.')   
       
    % Create display for the strings
    S_r = 20;
    y_p = 290;

    ah3 = axes('Parent',fh,'units','pixels',...
               'Position',[0 y_p res_size(1) 61 ],...
               'Visible', 'on',...
               'DataAspectRatioMode','auto'...
           ); 
    axis off
    set(ah3,'xlim',[0 800]);
    set(ah3,'ylim',[0 60]);

    [th,ws] = create_display_string(alignment,ah3,0, S_r, res_size);
    
    % Do slider for roling window
    SliderW = uicontrol('Parent', fh, 'Style', 'slider',...
                   'Position', [20 y_p-30 res_size(1)-40  20],...
                   'Min', 0, ...
                   'Max', ws - 800, ...
                   'Tag', 'SliderWindow', ...
                   'SliderStep', [0.05 0.2] ,...
                   'Callback',{@popup_menu_Callback} );   

    % show scores
    uicontrol('Style','text', ...
                    'String','Score:',...
                    'Position', [330 120 160 40],'FontSize', 18);
    Scorea = uicontrol('Style','text', ...
                    'String',score,...
                    'Position', [330 40 160 60], 'ForegroundColor', [1 0 0], 'FontSize', 24);



    % Display ScanMatchInfo parameters  
    ah1 = uipanel('units','pixels',...
               'Position',[520 0 279 239]...
           );

    uicontrol(ah1,'Style','text', ...
                    'String','Parameters',...
                    'Position', [0 210 270 25], 'FontSize', 16);

    uicontrol(ah1,'Style','text', ...
                    'String','Alignment algorithm:',...
                    'Position', [5 170 130 20], 'FontSize', 10, ...
                    'HorizontalAlignment', 'right');

    uicontrol(ah1,'Style','text', ...
                    'String','Needleman-Wunsch',...
                    'Position', [145 170 130 20], 'FontSize', 10, ...
                    'FontWeight', 'bold', 'HorizontalAlignment', 'left');

    uicontrol(ah1,'Style','text', ...
                    'String','Gap value:',...
                    'Position', [5 150 130 20], 'FontSize', 10, ...
                    'HorizontalAlignment', 'right');

    uicontrol(ah1,'Style','text', ...
                    'String',num2str(ScanMatchInfo.GapValue),...
                    'Position', [145 150 130 20], 'FontSize', 10, ...
                    'FontWeight', 'bold','HorizontalAlignment', 'left');

    uicontrol(ah1,'Style','text', ...
                    'String','Substitution threshold:',...
                    'Position', [5 130 130 20], 'FontSize', 10, ...
                    'HorizontalAlignment', 'right');

    uicontrol(ah1,'Style','text', ...
                    'String',num2str(ScanMatchInfo.Threshold),...
                    'Position', [145 130 130 20], 'FontSize', 10, ...
                    'FontWeight', 'bold','HorizontalAlignment', 'left');

    uicontrol(ah1,'Style','text', ...
                    'String','RoI modulus:',...
                    'Position', [5 110 130 20], 'FontSize', 10, ...
                    'HorizontalAlignment', 'right');

    uicontrol(ah1,'Style','text', ...
                    'String',num2str(ScanMatchInfo.RoiModulus),...
                    'Position', [145 110 130 20], 'FontSize', 10, ...
                    'FontWeight', 'bold','HorizontalAlignment', 'left');    

    uicontrol(ah1,'Style','text', ...
                    'String','Temporal Bin size:',...
                    'Position', [5 90 130 20], 'FontSize', 10, ...
                    'HorizontalAlignment', 'right');

    uicontrol(ah1,'Style','text', ...
                    'String',[num2str(ScanMatchInfo.TempBin) ' ms'],...
                    'Position', [145 90 130 20], 'FontSize', 10, ...
                    'FontWeight', 'bold','HorizontalAlignment', 'left');   

    uicontrol(ah1,'Style','text', ...
                    'String','Stimuli Resolution:',...
                    'Position', [5 70 130 20], 'FontSize', 10, ...
                    'HorizontalAlignment', 'right');

    uicontrol(ah1,'Style','text', ...
                    'String',[num2str(ScanMatchInfo.Xres) ' x ' num2str(ScanMatchInfo.Yres) ' pixels'],...
                    'Position', [145 70 130 20], 'FontSize', 10, ...
                    'FontWeight', 'bold','HorizontalAlignment', 'left');   

    uicontrol(ah1,'Style','text', ...
                    'String','Number of bins:',...
                    'Position', [5 50 130 20], 'FontSize', 10, ...
                    'HorizontalAlignment', 'right');

    uicontrol(ah1,'Style','text', ...
                    'String',[num2str(ScanMatchInfo.Xbin) ' x ' num2str(ScanMatchInfo.Ybin) ],...
                    'Position', [145 50 130 20], 'FontSize', 10, ...
                    'FontWeight', 'bold','HorizontalAlignment', 'left'); 
end
            
function popup_menu_Callback(source,eventdata) 
% This function is Nested hence all variables are globals. 
      pop_menu = get(source,'Tag');
      switch pop_menu
          case 'SliderGap'
              set(fh,'CurrentAxes',ah2)
              GapValue = get(source,'Value');
              [score, data]=nwalignds(StringA, StringB,'ScoringMatrix',Transition_dist, 'GapOpen', GapValue,'SHOWSCORE', true);
              set(Scorea, 'string', num2str(score));
              set(Gapa, 'string', ['Gap value: ' num2str(GapValue)]);
              set(fh,'CurrentAxes',ah3)
              delete(th); % delete old display
              [th ws] = create_display_string(data,ah3,0, S_r,info); % Create new display for the strings
              set(SliderW, 'Value', 0);
              set(SliderW, 'Max', ws-800);
          case 'Transition_mat'
              
          case 'SliderWindow'
              dx = get(source,'Value');
              set(ah3,'xlim',[dx 800+dx]);
              
      end
end

end
   

function [th, ws] = create_display_string(data,ah3,y_p, S_r, res_size)

n_string = size(data,2)/2; % IF MORE THAN 30 YOU SHOULD SHOW A SCROOLING BAR!
if n_string > 25 
   S_c = res_size(1)/25;
else
    S_c = res_size(1)/n_string-1;
end
ws = n_string * S_c;
th = zeros(3,n_string);
tic
for i=1:3
    x_p = S_c/2;
    ind=1;
    for j=1:n_string
    switch data(2,ind)
        case '|'
            col = [0 0 0];
        case ':'
            col = [1 0 0];
        case ' '
            col = [0 0 1];
    end
   
    th(i,j)=text('parent', ah3', 'Position', [x_p y_p], ...
            'string', data(i,ind:ind+1), 'color', col,'HorizontalAlignment', 'center' );
    x_p = x_p + S_c;
    ind = ind +2;
    end
    y_p = y_p + S_r;
end

toc
end

function pcmap = privateColorMap(selection)
%PRIVATECOLORMAP returns a custom color map
switch selection
    case 1, pts = [0 0 .3 20;
            0 .1 .8 25;
            0 .9 .5 15;
            .9 1 .9 8;
            1 1 0 26;
            1 0 0 26;
            .4 0 0 0];
    otherwise, pts = [0 0 0 128; 1 1 1 0];
end
xcl=1;
for i=1:size(pts,1)-1
    xcl=[xcl,i+1/pts(i,4):1/pts(i,4):i+1];
end
pcmap = interp1(pts(:,1:3),xcl);
end