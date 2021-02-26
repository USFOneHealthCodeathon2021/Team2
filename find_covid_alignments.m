function [reps] = find_covid_alignments(seq1,seq2)

% MAKE SURE TO ASSIGN SEQ1 AS WUHAN SEQUENCE AND REMOVE ARGUMENT

%~~~~~~~~~~~~~~~~~~~~~~~~

% optimum parameters for covid viruses:
min_len_inrep = 10;
match_len = 5;
indel_len = 0;

%~~~~~~~~~~~~~~~~~~~~~~~~

% find nonoverlapping alignments, with both sequences first for
% more completeness:
reps1 = make_align(seq1,seq2,min_len_inrep,match_len,indel_len);
reps2 = make_align(seq2,seq1,min_len_inrep,match_len,indel_len);

% combine the overlaps into one variables "reps":
reps = zeros(size(reps1,1)+size(reps2,1),4);
for r = 1:size(reps1,1)
    reps(r,1) = reps1{r,1}(1);
    reps(r,2) = reps1{r,1}(2);
    reps(r,3) = reps1{r,2}(1);
    reps(r,4) = reps1{r,2}(2);
end
for r = 1:size(reps2,1)
    reps(size(reps1,1)+r,1) = reps2{r,2}(1);
    reps(size(reps1,1)+r,2) = reps2{r,2}(2);
    reps(size(reps1,1)+r,3) = reps2{r,1}(1);
    reps(size(reps1,1)+r,4) = reps2{r,1}(2);
end

% merge the overlaps in reps:
r1 = 1;
while r1 <= size(reps,1)-1
    r2 = r1+1;
    while r2 <= size(reps,1)
        if reps(r1,1)==reps(r2,1) && reps(r1,2)==reps(r2,2) && reps(r1,3)==reps(r2,3) && reps(r1,4)==reps(r2,4)
            reps(r2,:) = [];
        elseif reps(r2,1) >= reps(r1,1) && reps(r2,1) <= reps(r1,2) || reps(r2,2) >= reps(r1,1) && reps(r2,2) <= reps(r1,2)
            if reps(r1,1)-reps(r2,1) == reps(r1,3)-reps(r2,3) && reps(r1,2)-reps(r2,2) == reps(r1,4)-reps(r2,4)
                reps(r1,1) = min(reps(r1,1),reps(r2,1));
                reps(r1,2) = max(reps(r1,2),reps(r2,2));
                reps(r1,3) = min(reps(r1,3),reps(r2,3));
                reps(r1,4) = min(reps(r1,4),reps(r2,4));
                reps(r2,:) = [];
            end
        else
            r2 = r2 + 1;
        end
    end
    r1 = r1 + 1;
end


% make a visual display of the sequences:
figure
hold on

title('Comparison of input CoVid sequence and original Wuhan strain sequence')
xlabel('base pairs')

YTick_vals = [0.5, 2.5];
YTick_labels = {'original Wuhan strain','input sequence'};
set(gca,'YTick',YTick_vals)
set(gca,'YTickLabels',YTick_labels)

len1 = length(seq1);
len2 = length(seq2);

spacer = 500;
loc = reps;
for r = 1:size(reps,1)
    loc(r,:) = loc(r,:) + spacer*r;
end

axis([0 max(loc(size(reps,1),2),loc(size(reps,1),4))+spacer 0 3])

for r = 1:size(reps,1)
    L1 = loc(r,1);
    R1 = loc(r,2);
    midpt1 = L1+(R1-L1)/2;
    L2 = loc(r,3);
    R2 = loc(r,4);
    midpt2 = L2+(R2-L2)/2;
    fill([R1,L1,L1,R1],[1,1,0,0],'g','FaceAlpha',0.2,'EdgeColor',[0 0.4 0])
    fill([R2,L2,L2,R2],[3,3,2,2],'g','FaceAlpha',0.2,'EdgeColor',[0 0.4 0])
    line([midpt1,midpt2],[1,2],'Color',[0 0.4 0],'LineStyle',':')
    text(midpt1,0.5,[num2str(R1-L1+1),' bp'],'HorizontalAlignment','center','Rotation',90)
    text(midpt2,2.5,[num2str(R2-L2+1),' bp'],'HorizontalAlignment','center','Rotation',90)
    
    if r == 1
        if reps(1,1) > 1
            misalign1 = seq1(1:reps(1,1)-1);
            if reps(1,1)-1 == 1
                bp = ['(1)'];
            else
                bp = ['(1:',num2str(reps(1,1)-1),')'];
            end
            text(spacer/2,0.5,[misalign1,' ',bp],'Rotation',90,'HorizontalAlignment','center')
        end
        if reps(1,3) > 1
            misalign2 = seq2(1:reps(1,3)-1);
            if reps(1,3)-1 == 1
                bp = '(1)';
            else
                bp = ['(1:',num2str(reps(1,3)-1),')'];
            end
            text(spacer/2,2.5,[misalign2,' ',bp],'Rotation',90,'HorizontalAlignment','center')
        end
    end
    if r < size(reps,1)
        misalign1 = seq1(reps(r,2)+1:reps(r+1,1)-1);
        if reps(r+1,1)-reps(r,2) == 2
            bp = ['(',num2str(reps(r,2)+1),')'];
        else
            bp = ['(',num2str(reps(r,2)+1),':',num2str(reps(r+1,1)-1),')'];
        end
        text(loc(r,2)+1+(loc(r+1,1)-1-(loc(r,2)+1))/2,0.5,[misalign1,' ',bp],'Rotation',90,'HorizontalAlignment','center')
        misalign2 = seq2(reps(r,4)+1:reps(r+1,3)-1);
        if reps(r+1,3)-reps(r,4) == 2
            bp = ['(',num2str(reps(r,4)+1),')'];
        else
            bp = ['(',num2str(reps(r,4)+1),':',num2str(reps(r+1,3)-1),')'];
        end
        text(loc(r,4)+1+(loc(r+1,3)-1-(loc(r,4)+1))/2,2.5,[misalign2,' ',bp],'Rotation',90,'HorizontalAlignment','center')
    elseif r == size(reps,1)
        if reps(size(reps,1),2) < len1
            misalign1 = seq1(reps(size(reps,1),2)+1:len1);
            if reps(r,2) == len1-1
                bp = ['(',num2str(len1),')'];
            else
                bp = ['(',num2str(reps(r,2)+1),':',num2str(len1),')'];
            end
            text(loc(size(loc,1),2)+spacer/2,0.5,[misalign1,' ',bp],'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom')
        end
        if reps(size(reps,1),2) < len1
            misalign2 = seq2(reps(size(reps,1),4)+1:len2);
            if reps(r,4) == len2-1
                bp = ['(',num2str(len2),')'];
            else
                bp = ['(',num2str(reps(r,4)+1),':',num2str(len2),')'];
            end
            text(loc(size(loc,1),4)+spacer/2,2.5,[misalign2,' ',bp],'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom')
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [reps] = make_align(seq1,seq2,min_len_inrep,match_len,indel_len)


% INPUTS: 
% seq1: a [1xm] character vector of nucleotides
% seq2: a [1xn] character vector of nucleotides

% OUTPUTS:
% reps: a [ix2] cell array, where i is the number of aligned sequences; 
% each of the i rows has format:
%   reps{i,1}(1:2): the range on seq1 of an aligned sequence
%   reps{i,2}(1:2): the range on seq2 of corresponding aligned sequence

% NOTE:
% for inverted repeats, the contents of inreps can be translated to the
% contig sequence in the following way, where the lower_range is the
% sequence to the 5' side of the area of interest, and upper_range is the
% sequence to the 3' side:
%   start of forward inverted repeat = lower_range(1)+inreps{i,1}(1)-1
%   end of forward inverted repeat = lower_range(1)+inreps{i,1}(2)-1
%   start of reverse inverted repeat = upper_range(2)-inreps{i,2}(2)+1
%   end of reverse inverted repeat = upper_range(2)-inreps{i,2}(1)+1

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

len1 = length(seq1);
len2 = length(seq2);

align1 = [];
align2 = [];

% Find the ranges for sequences of length match_len or longer that match:
OS1 = 0;
while OS1 <= len1 - match_len
    OS2 = 0;
    while OS2 <= len2 - match_len && OS1 <= len1 - match_len
        S1 = seq1(1+OS1:match_len+OS1);  %%%%%
        S2 = seq2(1+OS2:match_len+OS2);
        
        if strcmp(S1,S2)
            expanding = true; expand = 0;
            while expanding
                if OS1+match_len+expand+1 <= len1 && OS2+match_len+expand+1 <= len2 && strcmp(seq1(OS1+match_len+expand+1),seq2(OS2+match_len+expand+1))
                    expand = expand+1;
                else
                    expanding = false;
                end
            end
            % expand as much as possible
            align1 = [align1; 1+OS1, match_len+OS1+expand];
            align2 = [align2; 1+OS2, match_len+OS2+expand];
            
            % update OS1 and OS2 to be at the end of the aligned area:
            OS1 = match_len+OS1+expand+1;
            OS2 = match_len+OS2+expand+1;
        else
            OS2 = OS2 + 1;
        end
    end
    OS1 = OS1 + 1;
end


% Make a matrix of possible transitions between matching ranges:
aligns = size(align1,1);
mat = zeros(aligns);
for r = 1:aligns
    for c = 1:aligns
        diff1 = align1(c,1)-align1(r,2)-1;
        diff2 = align2(c,1)-align2(r,2)-1;
        if abs(diff1) <= indel_len && abs(diff2) <= indel_len && abs(diff1-diff2) <= indel_len
            mat(r,c) = 1;
        else
            mat(r,c) = 0;
        end
    end
end
 

% Find all possible starting ranges for longer sequences, including lone
% ranges:
starts = [];
for a = 1:aligns
    row_sum = sum(mat(a,:));
    col_sum = sum(mat(:,a));
    if row_sum > 0 && col_sum == 0 || row_sum == 0 && col_sum == 0
        starts = [starts,a];
    end
end


% Make a list of all possible sequences of ranges:
if ~isempty(starts)
    for s = 1:length(starts)
        listo(s).seq = starts(s);
    end
    k = 1;
    while k <= size(listo,2)
        len_seq = length(listo(k).seq);
        last_el = listo(k).seq(len_seq);
        tails = find(mat(last_el,:) ~= 0);
        tail = 1;
        while tail <= length(tails)
            if tails(tail) < last_el
                tails(tail) = [];
            else
                tail = tail + 1;
            end
        end
        if ~isempty(tails)
            listo(k).seq = [listo(k).seq, tails(1)];
            if length(tails) > 1
                for t = 2:length(tails)
                    listo(size(listo,2)+1).seq = [listo(k).seq(1:len_seq), tails(t)];
                end
            end
        else
            k = k + 1;
        end
    end
else
    listo = [];
end


% Delete sequences that are shorter than min_len_inrep and pack into
% variable inreps:
reps = {};
i = 0;
if ~isempty(listo)
    for k = 1:size(listo,2)
        first_el = listo(k).seq(1);
        last_el = listo(k).seq(size(listo(k).seq,2));
        range1 = [align1(first_el,1), align1(last_el,2)];
        range2 = [align2(first_el,1), align2(last_el,2)];
        
        if range1(2)-range1(1)+1 >= min_len_inrep && range2(2)-range2(1)+1 >= min_len_inrep
            i = i + 1;
            reps{i,1} = range1;
            reps{i,2} = range2;
        end 
    end
end

