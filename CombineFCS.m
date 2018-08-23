function CombineFCS(FCS1,FCS2, outputfilename)
% This function can be used to combine 2 FCS files having a set of shared
% markers and return one FCS file with the total number of cells is equal
% to the summation of cells in both FCS files, with each cell has an
% extended number of measured markers.
%
% Input description:
% FCS1: Name (full name with directory) of the first FCS file
%
% FCS2: Name (full name with directory) of the second FCS file
%
% outputfilename: It is the desired file name for the combined output FCS
%                 file, with the '.fcs' extension.
%
% Important Notes: 1) A small windows listing the variables will appear for
%                  each FCS inout file, in order to select only the markers
%                  used for analysis, excluding others like (Time_event,
%                  DNA, Sample_Tag, etc). It is not to select the shared
%                  markers set, the shared markers will be selected
%                  automatically.
%
%                  2) The shared markers short names (PnN or name2) must be 
%                  the same in both FCS files, in order to the function to
%                  identify them and use them for combination.
%
% For citation and further information please refer to this publication:
% "CyTOFmerge: CyTOFmerge: Integrating mass cytometry data across multiple panels"

% get the Markers names of both files
[fcsdat,fcshdr,~,~] = fca_readfcs(FCS1);
VarNames1=cellstr(char(fcshdr.par.name2));
checked1 = SelectionFigure(VarNames1, 1);
Data1=fcsdat(:,checked1);
VarNames1=VarNames1(checked1);

[fcsdat,fcshdr,~,~] = fca_readfcs(FCS2);
VarNames2=cellstr(char(fcshdr.par.name2));
checked2 = SelectionFigure(VarNames2, 2);
Data2=fcsdat(:,checked2);
VarNames2=VarNames2(checked2);

%find shared Markers
Matches=zeros(length(VarNames1),1);
for i=1:length(VarNames1)
    if (nnz(strcmp(VarNames2,VarNames1(i)))==0)
        Matches(i)=0;
    else
    Matches(i)=find(strcmp(VarNames2,VarNames1(i)));
    end
end

%Reorder Data
Shared_Index=find(Matches>0);
Data1_NonShared=Data1;
Data1_NonShared(:,Shared_Index)=[];
VarNames1_NonShared=VarNames1;
VarNames1_NonShared(Shared_Index)=[];
Data2_NonShared=Data2;
Data2_NonShared(:,Matches(Shared_Index))=[];
VarNames2_NonShared=VarNames2;
VarNames2_NonShared(Matches(Shared_Index))=[];

% This is the order data {Shared  Non_Shared}
Data1=[ Data1(:,Shared_Index) Data1_NonShared];
VarNames1= [VarNames1(Shared_Index); VarNames1_NonShared];
Data2=[ Data2(:,Matches(Shared_Index)) Data2_NonShared];
VarNames2= [VarNames2(Matches(Shared_Index)); VarNames2_NonShared];

% Combine Data
m = length(Shared_Index);      % Number of shared markers

[IDX1]=knnsearch(Data2(:,1:m),Data1(:,1:m),'K',50,'Distance','euclidean');
[IDX2]=knnsearch(Data1(:,1:m),Data2(:,1:m),'K',50,'Distance','euclidean');

Data_combine1 = zeros(size(Data1,1),size(Data1,2)+size(Data2,2)-m);
Data_combine2 = zeros(size(Data2,1),size(Data1,2)+size(Data2,2)-m);

for i=1:size(Data1,1)
    Data_combine1(i,1:m)=Data1(i,1:m);
    Data_combine1(i,m+1:size(Data1,2))=Data1(i,m+1:end);
    Data_combine1(i,size(Data1,2)+1:end)=median(Data2(IDX1(i,1:50),m+1:end),1);
end
clear i

for i=1:size(Data2,1)
    Data_combine2(i,1:m)=Data2(i,1:m);
    Data_combine2(i,m+1:size(Data1,2))=median(Data1(IDX2(i,1:50),m+1:end),1);
    Data_combine2(i,size(Data1,2)+1:end)=Data2(i,m+1:end);
end
clear i

Data_combine=[Data_combine1; Data_combine2];
VarNames_combine=[VarNames1; VarNames2(m+1:end)];
clear IDX1 IDX2 Data_combine1 Data_combine2

% create FCS header
header.BYTEORD='1,2,3,4';
header.DATATYPE='F';
header.MODE='L';
header.NEXTDATA=0;
header.PnN=VarNames_combine;
header.PnS=VarNames_combine;

writeFCS(outputfilename,Data_combine,header);
end
%%
function checked = SelectionFigure(VarNames, filenumber)

% Create figure
no_col = ceil(length(VarNames)/30);
h.f = figure('Name',['Select Markers List for file ' num2str(filenumber)],'NumberTitle','off',...
    'units','pixels','position',[300,50,150*no_col,670],...
    'toolbar','none','menu','none');

% Create checkboxes
for i=1:no_col
    for j=1:30
        if (((i-1)*30+j) <= length(VarNames))
            h.c((i-1)*30+j) = uicontrol('style','checkbox','units','pixels',...
                'position',[10+150*(i-1),640-(j-1)*20,140,20],'string',VarNames{(i-1)*30+j});
        end
    end
end
% Create OK pushbutton
h.p = uicontrol('style','pushbutton','units','pixels',...
    'position',[75*no_col-35,5,70,20],'string','OK',...
    'callback',@p_call);
checked = [];
waitfor(h.f);

% Pushbutton callback
    function p_call(varargin)
        vals = get(h.c,'Value');
        checked = find([vals{:}]);
        if isempty(checked)
            warndlg('Please select the measured markers','Warning');
        else
            close(h.f);
        end
    end

end
%%
function [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(filename)
% [fcsdat, fcshdr, fcsdatscaled] = fca_readfcs(filename);
%
% Read FCS 2.0 and FCS 3.0 type flow cytometry data file and put the list mode
% parameters to the fcsdat array with size of [NumOfPar TotalEvents].
% Some important header data are stored in the fcshdr structure:
% TotalEvents, NumOfPar, starttime, stoptime and specific info for parameters
% as name, range, bitdepth, logscale(yes-no) and number of decades.
%
% [fcsdat, fcshdr] = fca_readfcs;
% Without filename input the user can select the desired file
% using the standard open file dialog box.
%
% [fcsdat, fcshdr, fcsdatscaled] = fca_readfcs(filename);
% Supplying the third output the fcsdatscaled array contains the scaled
% parameters. It might be useful for logscaled parameters, but no effect
% in the case of linear parameters. The log scaling is the following
% operation for the "ith" parameter:
% fcsdatscaled(:,i) = ...
%   10.^(fcsdat(:,i)/fcshdr.par(i).range*fcshdr.par(i).decade;);
%
%
% Copyright (c) 2011, Laszlo Balkay
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
% Ver 2.5
% 2006-2009 / University of Debrecen, Institute of Nuclear Medicine
% Laszlo Balkay
% balkay@pet.dote.hu
%
% 14/08/2006 I made some changes in the code by the suggestion of
% Brian Harms <brianharms@hotmail.com> and Ivan Cao-Berg <icaoberg@cmu.edu>
% (given at the user reviews area of Mathwork File exchage) The program should work
% in the case of Becton EPics DLM FCS2.0, CyAn Summit FCS3.0 and FACSDiva type
% list mode files.
%
% 29/01/2008 Updated to read the BD LSR II file format and including the comments of
% Allan Moser (Cira Discovery Sciences, Inc.)
%
% 24/01/2009 Updated to read the Partec CyFlow format file. Thanks for
% Gavin A Price
%
% if noarg was supplied

if nargin == 0
    [FileName, FilePath] = uigetfile('*.*','Select fcs file');
    filename = [FilePath,FileName];
    if FileName == 0;
        fcsdat = []; fcshdr = [];
        return;
    end
else
    filecheck = dir(filename);
    if size(filecheck,1) == 0
        hm = msgbox([filename,': The file does not exist!'], ...
            'FcAnalysis info','warn');
        fcsdat = []; fcshdr = [];
        return;
    end
end

% if filename arg. only contain PATH, set the default dir to this
% before issuing the uigetfile command. This is an option for the "fca"
% tool
[FilePath, FileNameMain, fext] = fileparts(filename);
FilePath = [FilePath filesep];
FileName = [FileNameMain, fext];
if  isempty(FileNameMain)
    currend_dir = cd;
    cd(FilePath);
    [FileName, FilePath] = uigetfile('*.*','Select FCS file');
    filename = [FilePath,FileName];
    if FileName == 0;
        fcsdat = []; fcshdr = [];
        return;
    end
    cd(currend_dir);
end

%fid = fopen(filename,'r','ieee-be');
fid = fopen(filename,'r','b');
fcsheader_1stline   = fread(fid,64,'char');
fcsheader_type = char(fcsheader_1stline(1:6)');
%TMP: update to include FCS 3.1
if strcmp(fcsheader_type,'FCS3.1')
    fcsheader_type='FCS3.0';
end
%
%reading the header
%
if strcmp(fcsheader_type,'FCS1.0')
    hm = msgbox('FCS 1.0 file type is not supported!','FcAnalysis info','warn');
    fcsdat = []; fcshdr = [];
    fclose(fid);
    return;
elseif  strcmp(fcsheader_type,'FCS2.0') || strcmp(fcsheader_type,'FCS3.0') % FCS2.0 or FCS3.0 types
    fcshdr.fcstype = fcsheader_type;
    FcsHeaderStartPos   = str2num(char(fcsheader_1stline(11:18)'));
    FcsHeaderStopPos    = str2num(char(fcsheader_1stline(19:26)')); %RLF edited to full 8-byte length
    FcsDataStartPos     = str2num(char(fcsheader_1stline(27:34)')); %RLF edited to full 8-byte length
    status = fseek(fid,FcsHeaderStartPos,'bof');
    fcsheader_main = fread(fid,FcsHeaderStopPos-FcsHeaderStartPos+1,'char');%read the main header
    warning off MATLAB:nonIntegerTruncatedInConversionToChar;
    fcshdr.filename = FileName;
    fcshdr.filepath = FilePath;
    % "The first character of the primary TEXT segment contains the
    % delimiter" (FCS standard)
    if fcsheader_main(1) == 12
        mnemonic_separator = 'FF';
    elseif fcsheader_main(1) == 9
        mnemonic_separator = 'TAB'; %RLF
    else
        mnemonic_separator = char(fcsheader_main(1));
    end
    if mnemonic_separator == '@';% WinMDI
        hm = msgbox([FileName,': The file can not be read (Unsupported FCS type: WinMDI histogram file)'],'FcAnalysis info','warn');
        fcsdat = []; fcshdr = [];
        fclose(fid);
        return;
    end
    fcshdr.TotalEvents = str2num(get_mnemonic_value('$TOT',fcsheader_main, mnemonic_separator));
    if fcshdr.TotalEvents == 0
        fcsdat = 0;
        fcsdatscaled = 0;
        return
    end
    fcshdr.NumOfPar = str2num(get_mnemonic_value('$PAR',fcsheader_main, mnemonic_separator));
    fcshdr.Creator = get_mnemonic_value('CREATOR',fcsheader_main, mnemonic_separator);
    %comp matrix reader added by RLF 12_15_10
    comp = get_mnemonic_value('SPILLOVER',fcsheader_main,mnemonic_separator);
    if ~isempty(comp)
        %%%
        compcell=regexp(comp,',','split');
        nc=str2double(compcell{1});
        fcshdr.CompLabels=compcell(2:nc+1);
        fcshdr.CompMat=reshape(str2double(compcell(nc+2:end)'),[nc nc])';
    else
        fcshdr.CompLabels=[];
        fcshdr.CompMat=[];
    end
    plate = get_mnemonic_value('PLATE NAME',fcsheader_main,mnemonic_separator);
    if ~isempty(plate)
        fcshdr.plate=plate;
    end
    %%%%%%%%%%%%
    
    %%%%%%added by RLF to account for large files
    if FcsDataStartPos == 0
        FcsDataStartPos = str2num(get_mnemonic_value('$BEGINDATA',fcsheader_main, mnemonic_separator));
    end
    %%%%%%%%%%%%%%%%%%%%%
    
    for i=1:fcshdr.NumOfPar
        fcshdr.par(i).name = get_mnemonic_value(['$P',num2str(i),'N'],fcsheader_main, mnemonic_separator);
        fcshdr.par(i).name2 = get_mnemonic_value(['$P',num2str(i),'S'],fcsheader_main, mnemonic_separator);
        fcshdr.par(i).range = str2num(get_mnemonic_value(['$P',num2str(i),'R'],fcsheader_main, mnemonic_separator));
        fcshdr.par(i).bit = str2num(get_mnemonic_value(['$P',num2str(i),'B'],fcsheader_main, mnemonic_separator));
        %==============   Changed way that amplification type is treated ---  ARM  ==================
        par_exponent_str= (get_mnemonic_value(['$P',num2str(i),'E'],fcsheader_main, mnemonic_separator));
        if isempty(par_exponent_str)
            % There is no "$PiE" mnemonic in the Lysys format
            % in that case the PiDISPLAY mnem. shows the LOG or LIN definition
            islogpar = get_mnemonic_value(['P',num2str(i),'DISPLAY'],fcsheader_main, mnemonic_separator);
            if length(islogpar)==3 && isequal(islogpar, 'LOG')  % islogpar == 'LOG'
                par_exponent_str = '5,1';
            else % islogpar = LIN case
                par_exponent_str = '0,0';
            end
        end
        
        par_exponent= str2num(par_exponent_str);
        fcshdr.par(i).decade = par_exponent(1);
        if fcshdr.par(i).decade == 0
            fcshdr.par(i).log = 0;
            fcshdr.par(i).logzero = 0;
        else
            fcshdr.par(i).log = 1;
            if (par_exponent(2) == 0)
                fcshdr.par(i).logzero = 1;
            else
                fcshdr.par(i).logzero = par_exponent(2);
            end
        end
        gain_str = get_mnemonic_value(['$P',num2str(i),'G'],fcsheader_main, mnemonic_separator);
        if ~isempty(gain_str)
            fcshdr.par(i).gain=str2double(gain_str);
        else
            fcshdr.par(i).gain=1;
        end
        
        %============================================================================================
    end
    fcshdr.starttime = get_mnemonic_value('$BTIM',fcsheader_main, mnemonic_separator);
    fcshdr.stoptime = get_mnemonic_value('$ETIM',fcsheader_main, mnemonic_separator);
    fcshdr.cytometry = get_mnemonic_value('$CYT',fcsheader_main, mnemonic_separator);
    fcshdr.date = get_mnemonic_value('$DATE',fcsheader_main, mnemonic_separator);
    fcshdr.byteorder = get_mnemonic_value('$BYTEORD',fcsheader_main, mnemonic_separator);
    fcshdr.datatype = get_mnemonic_value('$DATATYPE',fcsheader_main, mnemonic_separator);
    fcshdr.system = get_mnemonic_value('$SYS',fcsheader_main, mnemonic_separator);
    fcshdr.project = get_mnemonic_value('$PROJ',fcsheader_main, mnemonic_separator);
    fcshdr.experiment = get_mnemonic_value('$EXP',fcsheader_main, mnemonic_separator);
    fcshdr.cells = get_mnemonic_value('$Cells',fcsheader_main, mnemonic_separator);
    fcshdr.creator = get_mnemonic_value('CREATOR',fcsheader_main, mnemonic_separator);
    fcshdr.cytsn = get_mnemonic_value('$CYTSN',fcsheader_main, mnemonic_separator);
else
    hm = msgbox([FileName,': The file can not be read (Unsupported FCS type)'],'FcAnalysis info','warn');
    fcsdat = []; fcshdr = [];
    fclose(fid);
    return;
end
%
%reading the events
%
status = fseek(fid,FcsDataStartPos,'bof');
if strcmp(fcsheader_type,'FCS2.0')
    if strcmp(mnemonic_separator,'\') || strcmp(mnemonic_separator,'FF')... %ordinary or FacsDIVA FCS2.0
            || strcmp(mnemonic_separator,'/') || strcmp(mnemonic_separator,'TAB')% added by GAP 1/22/09 %added by RLF 09/02/10
        if fcshdr.par(1).bit == 16
            fcsdat = uint16(fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16')');
            if strcmp(fcshdr.byteorder,'1,2')...% this is the Cytomics data
                    || strcmp(fcshdr.byteorder, '1,2,3,4') %added by GAP 1/22/09
                fcsdat = bitor(bitshift(fcsdat,-8),bitshift(fcsdat,8));
            end
        elseif fcshdr.par(1).bit == 32
            if fcshdr.datatype ~= 'F'
                fcsdat = (fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint32')');
            else % 'LYSYS' case
                fcsdat = (fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32')');
            end
        else
            bittype = ['ubit',num2str(fcshdr.par(1).bit)];
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],bittype, 'ieee-le')';
        end
    elseif strcmp(mnemonic_separator,'!');% Becton EPics DLM FCS2.0
        fcsdat_ = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16', 'ieee-le')';
        fcsdat = zeros(fcshdr.TotalEvents,fcshdr.NumOfPar);
        for i=1:fcshdr.NumOfPar
            bintmp = dec2bin(fcsdat_(:,i));
            fcsdat(:,i) = bin2dec(bintmp(:,7:16)); % only the first 10bit is valid for the parameter
        end
    end
    fclose(fid);
elseif strcmp(fcsheader_type,'FCS3.0')
    %     if strcmp(mnemonic_separator,'|') % CyAn Summit FCS3.0
    %         fcsdat_ = (fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16','ieee-le')');
    %         fcsdat = zeros(size(fcsdat_));
    %         new_xrange = 1024;
    %         for i=1:fcshdr.NumOfPar
    %             fcsdat(:,i) = fcsdat_(:,i)*new_xrange/fcshdr.par(i).range;
    %             fcshdr.par(i).range = new_xrange;
    %         end
    %     else % ordinary FCS 3.0
    %%%%%edited by RLF 06_30_10
    if strcmp(fcshdr.datatype,'D')
        if strcmp(fcshdr.byteorder, '1,2,3,4')
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'double','l')';
        elseif strcmp(fcshdr.byteorder,'4,3,2,1')
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'double','b')';
        end
    elseif strcmp(fcshdr.datatype,'F')
        if strcmp(fcshdr.byteorder, '1,2,3,4')
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32','l')';
        elseif strcmp(fcshdr.byteorder,'4,3,2,1')
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32','b')';
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %     end
    fclose(fid);
end
%
%calculate the scaled events (for log scales) %RLF added gain division
if nargout>2
    fcsdatscaled = zeros(size(fcsdat));
    for  i = 1 : fcshdr.NumOfPar
        Xlogdecade = fcshdr.par(i).decade;
        XChannelMax = fcshdr.par(i).range;
        Xlogvalatzero = fcshdr.par(i).logzero;
        if fcshdr.par(i).gain~=1
            fcsdatscaled(:,i)  = double(fcsdat(:,i))./fcshdr.par(i).gain;
            
        elseif fcshdr.par(i).log
            fcsdatscaled(:,i) = Xlogvalatzero*10.^(double(fcsdat(:,i))/XChannelMax*Xlogdecade);
        else fcsdatscaled(:,i)  = fcsdat(:,i);
        end
    end
    
end

if nargout>3 && ~isempty(fcshdr.CompLabels) %RLF. applied to fcsdatscaled rather than fcsdat.
    
    compcols=zeros(1,nc);
    colLabels={fcshdr.par.name};
    for i=1:nc
        compcols(i)=find(strcmp(fcshdr.CompLabels{i},colLabels));
    end
    fcsdatcomp=fcsdatscaled;
    fcsdatcomp(:,compcols)=fcsdatcomp(:,compcols)/fcshdr.CompMat;
else fcsdatcomp=[];
end
end
%%%%
function mneval = get_mnemonic_value(mnemonic_name,fcsheader,mnemonic_separator)

if strcmp(mnemonic_separator,'\')  || strcmp(mnemonic_separator,'!') ...
        || strcmp(mnemonic_separator,'|') || strcmp(mnemonic_separator,'@')...
        || strcmp(mnemonic_separator, '/') % added by GAP 1/22/08
    mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
    if isempty(mnemonic_startpos)
        mneval = [];
        return;
    end
    mnemonic_length = length(mnemonic_name);
    mnemonic_stoppos = mnemonic_startpos + mnemonic_length;
    next_slashes = findstr(char(fcsheader(mnemonic_stoppos+1:end)'),mnemonic_separator);
    next_slash = next_slashes(1) + mnemonic_stoppos;
    
    mneval = char(fcsheader(mnemonic_stoppos+1:next_slash-1)');
elseif strcmp(mnemonic_separator,'FF')
    mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
    if isempty(mnemonic_startpos)
        mneval = [];
        return;
    end
    mnemonic_length = length(mnemonic_name);
    mnemonic_stoppos = mnemonic_startpos + mnemonic_length ;
    next_formfeeds = find( fcsheader(mnemonic_stoppos+1:end) == 12);
    next_formfeed = next_formfeeds(1) + mnemonic_stoppos;
    
    mneval = char(fcsheader(mnemonic_stoppos + 1 : next_formfeed-1)');
elseif strcmp(mnemonic_separator,'TAB') %added by RLF August 2010
    mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
    if isempty(mnemonic_startpos)
        mneval = [];
        return;
    end
    mnemonic_length = length(mnemonic_name);
    mnemonic_stoppos = mnemonic_startpos + mnemonic_length ;
    next_formfeeds = find( fcsheader(mnemonic_stoppos+1:end) == 9);
    next_formfeed = next_formfeeds(1) + mnemonic_stoppos;
    
    mneval = char(fcsheader(mnemonic_stoppos + 1 : next_formfeed-1)');
end
end
%%
function TEXT = writeFCS(fname, DATA, TEXT, OTHER)
%writeFCS Save numeric data into FCS format used in flow cytometry.
%   writeFCS(FNAME, DATA) creates a ver 3.1 FCS file with file name FNAME
%       from the numerical matrix data. The program assings the parameters
%       by analyzing the properties of the numerical data.
%
%       FNAME               File name preferably with .fcs extension
%       DATA                NxM matrix with data. N matches the number of
%                           events (cells) and M represents the number of
%                           channels. All negative values will be converted
%                           to zero. The datatype of the matric will
%                           determine the datatype in the file. If data is
%                           integer only, it will be saved as integer,
%                           otherwise single or double will be used,
%                           depending on the datatype.
%
%   writeFCS(FNAME, DATA, TEXT)  creates an FCS file with file name FNAME
%       from the numerical matrix data. The program takes parameters from
%       the struct TEXT or if they are missing, it guesses them from the
%       properties of the numerical data.
%
%       TEXT                is a struct with fields defined by the "Data
%                           File Standard for Flow Cytometry Version FCS
%                           3.1 Normative Reference". The fields can
%                           contain either text or numbers depending on
%                           their nature. If not all required fields are
%                           supplied, they will be guessed from the
%                           properties of the data.
%
%   writeFCS(FNAME, DATA, TEXT, OTHER)  creates an FCS file with file name
%       FNAME from the numerical matrix data. The program takes parameters
%       from the struct TEXT or if they are missing, it guesses them from
%       the properties of the numerical data. OTHER is a struct that can
%       append any other custom data to the end of the FCS file
%
%       OTHER               is a struct with fields defined by the "Data
%                           File Standard for Flow Cytometry Version FCS
%                           3.1 Normative Reference". The fields can be
%                           chosen arbitrarily as long as they adhere to
%                           the format specification.
%   Examples:
%
%       % This example creates an integer FCS file from random data points
%       DATA = round(horzcat(randn(1e4, 1)*20+120, randn(1e4, 1)*80+350));
%       writeFCS('integerFCS.fcs', DATA)
%
%       % This example creates an integer FCS file from random data points
%       DATA = single(horzcat(randn(1e4, 1)*20+120, randn(1e4, 1)*80+350));
%       writeFCS('singleFCS.fcs', DATA)
%
%
% Copyright 2013 Jakub Nedbal
% $Revision: 1.0 $  $Date: 2013/07/11 16:07:00 $

% check if TEXT is defined, else create a struct
if nargin < 3
    TEXT.BYTEORD = '1,2,3,4';
end

% If TEXT hasn't got BYTEORD specified, add it
if ~isfield(TEXT, 'BYTEORD')
    TEXT.BYTEORD = '1,2,3,4';
end

%% Parse the text, correct any errors if encountered and fix the DATA if
%  necessary according to the values in TEXT.
[DATA, TEXT] = parseTEXT(DATA, TEXT);


%% FCS standard version
fcsver = 3.0;


%% Convert filednames in TEXT into a consistent string
STR = fields2string(TEXT);

%% data offsets
Stext = 256;                            % Position of text start
Etext = Stext + numel(STR) - 1;         % Position of text end
Sdata = 2 ^ ceil(log2(Etext + 1));      % Position of data start
Edata = Sdata + size(DATA, 1) * sum(TEXT.PnB(:)) / 8 - 1;
% Position of data end
Sgate = 0; Egate = 0;                   % Position of gates
if nargin > 3
    Sother = 256 * ceil((Edata + 1) / 256);
    % Position of other start
else
    Sother = [];                        % Leave empty if not in use
end

%% create text
outText = sprintf('FCS%3.01f    %8d%8d%8d%8d%8d%8d%8d', ...
    fcsver, ...
    Stext, ...
    Etext, ...
    Sdata, ...
    Edata, ...
    Sgate, ...
    Egate, ...
    Sother);
outText = sprintf('%s%s%s', outText, repmat(' ', 1, Stext - numel(outText)), STR);
outText = sprintf('%s%s', outText, repmat(' ', 1, Sdata - Etext - 1));

%% conform to the endianness of the output data
fields = fieldnames(TEXT);
if strcmp(TEXT.(fields{cellfun(@(x) strcmpi(x, 'byteord'), fields)}), ...
        '4,3,2,1')
    endian = true;
    machineformat = 'b';
else
    endian = false;
    machineformat = 'l';
end

%% Create data to match the data type and bites per number
datatype = fields{cellfun(@(x) strcmpi(x, 'datatype'), fields)};
switch TEXT.(datatype)
    case 'F'
        DT = 'single';
        DATAtmp = single(DATA');
    case 'D'
        DT = 'double';
        DATAtmp = double(DATA');
    case 'I'
        DT = 'uint8';
        % The DATA needs to be separated into 8-bit numbers.
        DATAtmp = zeros(sum(TEXT.PnB(:)) / 8, size(DATA, 1), 'uint8');
        u = 0;
        for i = 1 : size(DATA, 2)
            tmp = round(double(DATA(:, i)));
            % conform to endianness of data
            in = 1 : TEXT.PnB(i) / 8;       % big endian
            if endian
                in = fliplr(in);            % little endian
            end
            for j = 1 : numel(in)
                DATAtmp(u + in(j), :) = mod(floor(tmp / (256 ^ (j - 1))), 256);
            end
            u = u + j;
        end
        
end

fid = fopen(fname, 'w');
fprintf(fid, outText);
%ftell(fid)

fwrite(fid, DATAtmp, DT, 0, machineformat);
%fwrite(fid, DATA', 'uint16');
% if OTHER segment of the FCS file is specified
if nargin > 3
    outText = repmat(' ', 1, Sother - Edata - 1);
    STR = fields2string(OTHER);
    outText = sprintf('%s%s', outText, STR);
    fprintf(fid, outText);
end

fprintf(fid, '\n');
%ftell(fid)
fclose(fid);

end     % end of writeFCS function

function STR = fields2string(TEXT)
%% Convert fieldnames in TEXT into a consistent string
dlm = '/';      % delimiter

fields = fieldnames(TEXT);
STR = dlm;

for i = 1 : numel(fields)
    value = TEXT.(fields{i});
    field = upper(fields{i});
    if isnumeric(value) || iscell(value)
        for j = 1 : numel(value)
            f = field;
            if regexp(f, '^((PN)|(GN)|(RN))')
                f = [f(1), num2str(j), f(3 : end)];
            end
            if isnumeric(value)
                v = num2str(value(j));
            else
                v = value{j};
            end
            STR = sprintf('%s$%s%s%s%s', STR, f, dlm, v, dlm);
        end
    else
        STR = sprintf('%s$%s%s%s%s', STR, field, dlm, value, dlm);
    end
end

end     % End of fields2string function

function [DATA, TEXT] = parseTEXT(DATA, TEXT)

%% check if data is numeric
if ~isnumeric(DATA)
    error('DATA must be a numeric matrix');
end

%% check if text is a struct
if ~isstruct(TEXT)
    warning('TEXT must be a struct. It will be recreated automatically.')
    clear TEXT;
    TEXT.BYTEORD = '1,2,3,4';
end

%% Get the fieldnames of TEXT and process them
fns = fieldnames(TEXT);

% Process each required fieldname

%% $BYTEORD Byte order for data acquisition computer.
al = {'1,2,3,4', '4,3,2,1'};
fn = fns{cellfun(@(x) strcmpi(x, 'BYTEORD'), fns)};
if isempty(fn)
    TEXT.(fn) = '1,2,3,4';
    fprintf('Setting $BYTEORD to "1,2,3,4".\n');
else
    va = TEXT.(fn);
    if ~ischar(va)
        va = 'mistake';
    end
    if ~any(strcmpi(va, al))
        TEXT.(fn) = '1,2,3,4';
        warning('$BYTEORD must be either "1,2,3,4" or "4,3,2,1"');
    end
end


%% $DATATYPE Type of data in DATA segment (ASCII, integer, floating point).
al = {'I', 'F', 'D'}; nal = {'A'};
in = cellfun(@(x) strcmpi(x, 'DATATYPE'), fns);
if ~in
    fn = 'DATATYPE';
    TEXT = finddatatype(DATA, TEXT, fn);
else
    fn = fns{in};
    va = TEXT.(fn);
    if ~ischar(va)
        va = 'mistake';
    end
    if ~any(strcmpi(va, al))
        warning(['$DATATYPE must be "I", "F" or "D"\nfor unsigned', ...
            'integer, single or double float, respectively.']);
        if ~any(strcmpi(va, nal))
            warning('writeFCS does not support $DATATYPE "A".\n');
        end
        [DATA, TEXT] = finddatatype(DATA, TEXT, fn);
    end
end
TEXT.(fn) = upper(TEXT.(fn));


%% $MODE Data mode (list mode - preferred, histogram - deprecated).
al = {'L'}; nal = {'C', 'U'};
in = cellfun(@(x) strcmpi(x, 'MODE'), fns);
if ~in
    fn = 'MODE';
    TEXT.(fn) = 'L';
    fprintf('Setting $MODE to "L".\n');
else
    fn = fns{in};
    va = TEXT.(fn);
    if ~ischar(va)
        va = 'mistake';
    end
    if ~any(strcmpi(va, al))
        warning('$MODE must be "L" for list mode.');
        if ~any(strcmpi(va, nal))
            warning('writeFCS does not support $MODE "C" or "U".\n');
        end
        TEXT.(fn) = 'L';
    end
end
TEXT.(fn) = upper(TEXT.(fn));


%% $NEXTDATA Byte offset to next data set in the file.
in = cellfun(@(x) strcmpi(x, 'NEXTDATA'), fns);
if ~in
    fn = 'NEXTDATA';
    TEXT.(fn) = 0;
    fprintf('Setting $NEXTDATA to 0.\n');
else
    fn = fns{in};
    va = TEXT.(fn);
    if ~isnumeric(va)
        TEXT.(fn) = 0;
        warning('writeFCS only supports $NEXTDATA equal to 0 (zero).');
    elseif va ~= 0
        TEXT.(fn) = 0;
        warning('writeFCS only supports $NEXTDATA equal to 0 (zero).');
    end
end


%% $PAR Number of parameters in an event.
in = cellfun(@(x) strcmpi(x, 'PAR'), fns);
if ~in
    fn = 'PAR';
    TEXT.(fn) = size(DATA, 2);
    fprintf('Setting $PAR to %d.\n', size(DATA, 2));
else
    fn = fns{in};
    va = TEXT.(fn);
    if ~isnumeric(va)
        TEXT.(fn) = size(DATA, 2);
        warning('writeFCS only supports numerical integer $PAR.');
        fprintf('Setting $PAR to %d.\n', size(DATA, 2));
    elseif va ~= size(DATA, 2);
        TEXT.(fn) = size(DATA, 2);
        warning('$PAR must be equal to the second dimension of DATA.');
        fprintf('Setting $PAR to %d.\n', size(DATA, 2));
    end
end
PAR = TEXT.(fn);


%% $TOT Total number of events in the data set.
in = cellfun(@(x) strcmpi(x, 'TOT'), fns);
if ~in
    fn = 'TOT';
    TEXT.(fn) = size(DATA, 1);
    fprintf('Setting $TOT to %d.\n', size(DATA, 1));
else
    fn = fns{in};
    va = TEXT.(fn);
    if ~isnumeric(va)
        TEXT.(fn) = size(DATA, 1);
        warning('writeFCS only supports numerical integer $TOT.');
        fprintf('Setting $TOT to %d.\n', size(DATA, 1));
    elseif va ~= size(DATA, 1);
        TEXT.(fn) = size(DATA, 1);
        warning('$TOT must be equal to the first dimension of DATA.');
        fprintf('Setting $TOT to %d.\n', size(DATA, 1));
    end
end


%% $PnB Number of bits reserved for parameter number n.
in = cellfun(@(x) strcmpi(x, 'PnB'), fns);
if in
    fn = fns{in};
    va = TEXT.(fn);
    if ~isnumeric(va)
        TEXT.(fn) = [];
        warning('writeFCS only supports numerical integer $PnB.');
    elseif numel(va) ~= PAR;
        TEXT.(fn) = [];
        warning(['Size of $PnB must be equal to the second ', ...
            'dimension of DATA.']);
    end
end
TEXT = checkPnB(DATA, TEXT);


%% $PnE Amplification type for parameter n.
in = cellfun(@(x) strcmpi(x, 'PnE'), fns);
if in
    fn = fns{in};
    va = TEXT.(fn);
    if ~(iscell(va) || (ischar(va) && PAR == 1))
        TEXT.(fn) = {};
        warning('writeFCS only supports cell or character $PnE.');
    end
end
TEXT = checkPnE(DATA, TEXT);


%% $PnN name for parameter n.
in = cellfun(@(x) strcmpi(x, 'PnN'), fns);
if ~in
    fn = 'PnN';
    va = cellfun(@(x) sprintf('Ch %d', x), num2cell(1 : size(DATA, 2)), ...
        'UniformOutput', false);
    err = true;
else
    fn = fns{in};
    va = TEXT.(fn);
    err = false;
end

% check if TEXT.PnN is defined correctly
if isempty(va) || (ischar(va) && size(DATA, 2) > 1) || ...
        (iscell(va) && size(DATA, 2) ~= numel(va))
    warning(['$PnE must be a cell. Its size must be equal to the ', ...
        'second dimension of DATA.']);
    va = cellfun(@(x) sprintf('Ch %d', x), num2cell(1 : size(DATA, 2)), ...
        'UniformOutput', false);
    err = true;
end

% if TEXT.PnN is character, convert to a cell
if ischar(va)
    va = {va};
end

% Check that all values are a character, replace those which aren't
in = find(~cellfun(@ischar, va));
va(in) = cellfun(@(x) sprintf('Ch %d', x), num2cell(in), ...
    'UniformOutput', false);
if ~isempty(in)
    warning(['$PnE must be a cell. Its size must be equal to the ', ...
        'second dimension of DATA.']);
    err = true;
end

% Make sure that commas are not used in the names, replace them with dashes
if any(~cellfun(@isempty, strfind(va, ',')))
    warning(['$PnE values must not include commas (,). Replacing by ', ...
        'dashes (-).']);
    err = true;
    va = cellfun(@(x) regexprep(x, ',', '-'), va, 'UniformOutput', false);
end

if any(err)
    txt = '';
    for i = 1 : numel(va)
        txt = sprintf('%s''%s'', ', txt, va{i});
    end
    fprintf('Setting $PnN to {%s}.\n', txt(1 : end - 2));
end
TEXT.(fn) = va;


%% $PnS Short name for parameter n.
in = cellfun(@(x) strcmpi(x, 'PnS'), fns);
if ~in
    fn = 'PnS';
    va = cellfun(@(x) sprintf('Ch %d', x), num2cell(1 : size(DATA, 2)), ...
        'UniformOutput', false);
    err = true;
else
    fn = fns{in};
    va = TEXT.(fn);
    err = false;
end

% check if TEXT.PnN is defined correctly
if isempty(va) || (ischar(va) && size(DATA, 2) > 1) || ...
        (iscell(va) && size(DATA, 2) ~= numel(va))
    warning(['$PnE must be a cell. Its size must be equal to the ', ...
        'second dimension of DATA.']);
    va = cellfun(@(x) sprintf('Ch %d', x), num2cell(1 : size(DATA, 2)), ...
        'UniformOutput', false);
    err = true;
end

% if TEXT.PnN is character, convert to a cell
if ischar(va)
    va = {va};
end

% Check that all values are a character, replace those which aren't
in = find(~cellfun(@ischar, va));
va(in) = cellfun(@(x) sprintf('Ch %d', x), num2cell(in), ...
    'UniformOutput', false);
if ~isempty(in)
    warning(['$PnE must be a cell. Its size must be equal to the ', ...
        'second dimension of DATA.']);
    err = true;
end

% Make sure that commas are not used in the names, replace them with dashes
if any(~cellfun(@isempty, strfind(va, ',')))
    warning(['$PnE values must not include commas (,). Replacing by ', ...
        'dashes (-).']);
    err = true;
    va = cellfun(@(x) regexprep(x, ',', '-'), va, 'UniformOutput', false);
end

if any(err)
    txt = '';
    for i = 1 : numel(va)
        txt = sprintf('%s''%s'', ', txt, va{i});
    end
    fprintf('Setting $PnS to {%s}.\n', txt(1 : end - 2));
end
TEXT.(fn) = va;


%% $PnR Range for parameter number n.
in = cellfun(@(x) strcmpi(x, 'PnR'), fns);
if ~in
    fn = 'PnR';
    va = 2 .^ ceil(log2(max(DATA)));
    err = true;
else
    fn = fns{in};
    va = TEXT.(fn);
    err = false;
end

% check if TEXT.PnR is defined correctly
if isempty(va) || ~isnumeric(va)
    if ~isnumeric(va)
        warning(['$PnR must be a numeric array. Its size must be ', ...
            'equal to the second dimension of DATA.']);
    end
    va = 2 .^ ceil(log2(max(DATA)));
    err = true;
end

% check if values of TEXT.PnR have integer values
in = va ~= round(va);
if any(in)
    err = true;
    va(in) = 2 .^ ceil(log2(max(DATA(in))));
    warning('$PnR must consist of integer values.');
end

% check if values of TEXT.PnR have positive values
in = va < 1;
if any(in)
    err = true;
    va(in) = 2 .^ ceil(log2(max(DATA(in))));
    warning('$PnR must consist of positive values.');
end


if err
    txt = sprintf('%d, ', va);
    fprintf('Setting $PnR to [%s].\n', txt(1 : end - 2));
end
TEXT.(fn) = va;

end     % end of parseTEXT function

function TEXT = finddatatype(DATA, TEXT, fn)
%% Work out the data type based on the properties of the DATA
% The functions checks whether data is integer, single or double. If single
% or double and consisting only of integer values, it converts to integer.

if isinteger(DATA)
    TEXT.(fn) = 'I';
    fprintf('Setting $DATATYPE to "I".\n');
elseif isfloat(DATA)
    % Check if all values are integers
    if ~any(DATA(:) - round(DATA(:)))
        TEXT.(fn) = 'I';
        fprintf('Setting $DATATYPE to "I".\n');
    else
        switch class(DATA)
            case 'single'
                TEXT.(fn) = 'F';
                fprintf('Setting $DATATYPE to "F".\n');
            case 'double'
                TEXT.(fn) = 'D';
                fprintf('Setting $DATATYPE to "D".\n');
            otherwise
                TEXT.(fn) = 'I';
                warning('Cannot determine $DATATYPE. Setting to "I".');
        end
    end
end

end     % end of finddatatype function

function TEXT = checkPnB(DATA, TEXT)
%% Check the values fo $PnB and make sure they comply with FCS format
%  specification.
fns = fieldnames(TEXT);
in = cellfun(@(x) strcmpi(x, 'PnB'), fns);
if ~in
    fn = 'PnB';
    va = [];
    TEXT.(fn) = [];
    err = false;
else
    fn = fns{in};
    va = TEXT.(fn);
    err = true;
end

switch TEXT.(fns{cellfun(@(x) strcmpi(x, 'DATATYPE'), fns)})
    case 'I'        % integer
        uints = [8; 16; 32; 64];
        in = ~any(repmat(va(:)', numel(uints), 1) == ...
            repmat(uints, 1, numel(va)));
        if isempty(va)
            in = true(1, size(DATA, 2));
        elseif any(in) && err
            warning(['With integer $DATATYPE, $PnB can only be ', ...
                '8, 16, 32 or 64']);
        end
        va(in) = 2 .^ sum(repmat(max(DATA(:, in)), 3, 1) > ...
            repmat(2 .^ [8; 16; 32], 1, sum(in))) * 8;
        
    case 'F'        % single float
        in = ~any(va(:)' == 32 * ones(1, numel(va)));
        if (any(in) || isempty(in)) && err
            warning('With single float $DATATYPE, $PnB can only be 32.');
        end
        va = 32 * ones(1, size(DATA, 2));
        
    case 'D'        % double float
        in = ~any(va(:)' == 64 * ones(1, numel(va)));
        if (any(in) || isempty(in)) && err
            warning('With single float $DATATYPE, $PnB can only be 32.');
        end
        va = 64 * ones(1, size(DATA, 2));
end

if ~isequal(va, TEXT.(fn))
    txt = sprintf('%d, ', va);
    fprintf('Setting $PnB to [%s].\n', txt(1 : end - 2));
    TEXT.(fn) = va;
end

end     % end of checkPnB function

function TEXT = checkPnE(DATA, TEXT)
fns = fieldnames(TEXT);
in = cellfun(@(x) strcmpi(x, 'PnE'), fns);
if ~in
    fn = 'PnE';
    va = repmat({''}, 1, TEXT.(fns{cellfun(@(x) strcmpi(x, 'PAR'), fns)}));
    TEXT.(fn) = va;
    % error vector (zeros mean no error), last bit =0 means no warnings
    err = false(1, 5);
else
    fn = fns{in};
    va = TEXT.(fn);
    % error vector (zeros mean no error), last bit =1 means allow warnings
    err = [false(1, 4), true];
end

% check if it is a characted array with a single-dimensional data array
if ischar(va) && size(DATA, 2) == 1
    va = {va};
end

% check that it is a cell, else redefine it
if ~iscell(va)
    va = repmat({''}, 1, TEXT.(fns{cellfun(@(x) strcmpi(x, 'PAR'), fns)}));
    TEXT.(fn) = va;
    warning('$PnE must be a cell.')
end

switch TEXT.(fns{cellfun(@(x) strcmpi(x, 'DATATYPE'), fns)})
    case 'I'        % integer
        for i = 1 : size(DATA, 2)
            if numel(va) >= i
                if isempty(va{i})
                    cva = {[], []};
                else
                    cva = textscan(va{i}, '%f,%f');
                end
                if any(cellfun(@isempty, cva))
                    cva = {0, 0};
                    err(1) = true;
                end
            else
                cva = {0, 0};
                err(2) = true;
            end
            if any(cellfun(@(x) x < 0, cva))
                cva = {0, 0};
                err(3) = true;
            end
            if cva{1} > 0 && cva{2} == 0
                cva{2} = 10 .^ floor(log10(min(DATA(:, i))));
                err(3) = true;
            end
            va{i} = sprintf('%g,%g', cell2mat(cva));
        end
        if all(err([1 end]))
            warning(['$PnE must consist of two numbers separated by ', ...
                'a comma such as "0,0" or "4.0,0.1".']);
        end
        if all(err([2 end]))
            warning(['Size of $PnE must be equal to the second ', ...
                'dimension of DATA.']);
        end
        if all(err([3 end]))
            warning('$PnE must consist of non-negative numbers only.');
        end
        if all(err([4 end]))
            warning(['Unless being "0,0", $PnE must consist of ', ...
                'positive numbers only.']);
        end
        
    case {'F', 'D'}     % single or double float
        if numel(va) ~= size(DATA, 2)
            err = true;
            warning(['Size of $PnE must be equal to the second ', ...
                'dimension of DATA.']);
        end
        if ~iscell(va)
            warning('$PnE must be a cell');
            err = true;
        end
        try
            if any(~cellfun(@(x) isequal({0, 0}, textscan(x, '%f,%f')), va))
                err = true;
                warning(['With float $DATATYPE "F" or "D", $PnE must ', ...
                    'remain "0,0".']);
            end
        catch
            err(1) = true;
        end
        va = repmat({'0,0'}, 1, size(DATA, 2));
end

if any(err(1 : end - 1))
    txt = '';
    for i = 1 : numel(va)
        txt = sprintf('%s''%s'', ', txt, va{i});
    end
    fprintf('Setting $PnE to {%s}.\n', txt(1 : end - 2));
end
TEXT.(fn) = va;

end     % end of checkPnE function
