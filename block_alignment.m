function [align1, align2, overlap1, overlap2, frac_overlap, frac_align] = block_alignment(seq1,seq2,NW_params,alignment_params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS COPY FOR CODE-A-THON %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


show_flags = true;

min_len_seqs = alignment_params(1);
match_len = alignment_params(2);
indel_len = alignment_params(3);

% find matching bits:
if show_flags
    disp('Finding alignments')
end
reps = find_nonoverlapping_alignments(seq1,seq2,min_len_seqs,match_len,indel_len);
num_aligns = size(reps,1);


% put reps into a matrix:
aligns = zeros(num_aligns,4);
for r = 1:num_aligns
    aligns(r,1) = reps{r,1}(1);
    aligns(r,2) = reps{r,1}(2);
    aligns(r,3) = reps{r,2}(1);
    aligns(r,4) = reps{r,2}(2);
end


% Find the best set of nonoverlapping alignments:
% Make a matrix of possible transitions between matching ranges:
if show_flags
    disp('Making a matrix of alignments')
end
num_aligns = size(aligns,1);
mat = zeros(num_aligns);
for r = 1:num_aligns
    for c = 1:num_aligns
        if aligns(r,2) < aligns(c,1) && aligns(r,4) < aligns(c,3)
            mat(r,c) = 1;
        else
            mat(r,c) = 0;
        end
    end
end
 

% Find all possible starting ranges for longer sequences, including lone
% ranges:
starts = [];
for a = 1:num_aligns
    row_sum = sum(mat(a,:));
    col_sum = sum(mat(:,a));
    if row_sum > 0 && col_sum == 0 || row_sum == 0 && col_sum == 0
        starts = [starts,a];
    end
end


% Make a list of all possible sequences of ranges:
if show_flags
    disp('Making a list of possible sets of ranges')
end
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


% Find the list that includes the most base pairs:
if ~isempty(listo)
    max_len = 0;
    for j = 1:size(listo,2)
        len = 0;
        for k = 1:size(listo(j).seq,2)
            len = len + aligns(listo(j).seq(k),2)-aligns(listo(j).seq(k),1)+1;
        end
        if len > max_len
            max_len = len;
            max_listo = listo(j).seq;
        end
    end
end


% replace the ranges in seq1 and seq2 that correspond to alignments:
if show_flags
    disp('Replacing sequence ranges with blocks')
end
unique_els = unique([seq1,seq2]);
num_unique_els = length(unique_els);
num_aligns = length(max_listo);

newseq1 = cellstr(seq1')';
for r = 1:size(newseq1,2)
    found = false;
    el = 1;
    while ~found
        if strcmp(newseq1{r},unique_els(el))
            found = true;
            newseq1{r} = el;
        else
            el = el + 1;
        end
    end
end
newseq1 = cell2mat(newseq1);

for r = num_aligns:-1:1
    newseq1(aligns(max_listo(r),1)) = r+num_unique_els;
    for s = aligns(max_listo(r),1)+1:aligns(max_listo(r),2)
        newseq1(aligns(max_listo(r),1)+1) = [];
    end
end

newseq2 = cellstr(seq2')';
for r = 1:size(newseq2,2)
    found = false;
    el = 1;
    while ~found
        if strcmp(newseq2{r},unique_els(el))
            found = true;
            newseq2{r} = el;
        else
            el = el + 1;
        end
    end
end
newseq2 = cell2mat(newseq2);
for r = num_aligns:-1:1
    newseq2(aligns(max_listo(r),3)) = r+num_unique_els;
    for s = aligns(max_listo(r),3)+1:aligns(max_listo(r),4)
        newseq2(aligns(max_listo(r),3)+1) = [];
    end
end


% Make a vector of match_val weights:
num_weights = size(max_listo,2);
match_vals = zeros(2,num_weights);
match_vals(1,:) = linspace(1,num_weights,num_weights)+num_unique_els;
for k = 1:num_weights
    match_vals(2,k) = aligns(max_listo(k),2)-aligns(max_listo(k),1)+1;
end


% Align using modified Needleman-Wunsch aligner (subfunction):
if show_flags
    red1 = (length(seq1)-length(newseq1))/length(seq1)*100;
    red2 = (length(seq2)-length(newseq2))/length(seq2)*100;
    announcement = ['seq1 reduced in length by ',num2str(red1),'%, and seq2 by ',num2str(red2),'%'];
    disp(announcement)
    disp('Using Needleman-Wunsch aligner')
end
[block_align1,block_align2] = NW_block_align(newseq1,newseq2,match_vals,num_unique_els,NW_params);


% Replace aligned contents with nucleotide bases and underscores:
if show_flags
    disp('Replacing blocks with nucleotides')
end


disp_blocks_1 = '';
disp_blocks_2 = '';
align1 = '';
align2 = '';
for a = 1:size(block_align1,2)
    if block_align1(a) <= num_unique_els && block_align2(a) <= num_unique_els
        if block_align1(a) == 0
            align1 = strcat(align1,'_');
            disp_blocks_1 = strcat(disp_blocks_1,'_');
        else
            align1 = strcat(align1,unique_els(block_align1(a)));
            disp_blocks_1 = strcat(disp_blocks_1,'_');
        end
        if block_align2(a) == 0
            align2 = strcat(align2,'_');
            disp_blocks_2 = strcat(disp_blocks_2,'_');
        else
            align2 = strcat(align2,unique_els(block_align2(a)));
            disp_blocks_2 = strcat(disp_blocks_2,'_');
        end
    elseif block_align1(a) > num_unique_els && block_align2(a) == 0
        align1 = strcat(align1,seq1(aligns(block_align1(a)-num_unique_els,1):aligns(block_align1(a)-num_unique_els,2)));
        num = aligns(block_align1(a)-num_unique_els,2)-aligns(block_align1(a)-num_unique_els,1)+1;
        align2 = strcat(align2,repmat('_',1,num));
        disp_blocks_1 = strcat(disp_blocks_1,repmat('|',1,num));
        disp_blocks_2 = strcat(disp_blocks_1,repmat('_',1,num));
    elseif block_align1(a) == 0 && block_align2(a) > num_unique_els
        align2 = strcat(align2,seq2(aligns(block_align2(a)-num_unique_els,3):aligns(block_align2(a)-num_unique_els,4)));
        num = aligns(block_align2(a)-num_unique_els,4)-aligns(block_align2(a)-num_unique_els,3)+1;
        align1 = strcat(align1,repmat('_',1,num));
        disp_blocks_1 = strcat(disp_blocks_1,repmat('_',1,num));
        disp_blocks_2 = strcat(disp_blocks_1,repmat('|',1,num));
    else
        block1 = seq1(aligns(block_align1(a)-num_unique_els,1):aligns(block_align1(a)-num_unique_els,2));
        block2 = seq2(aligns(block_align2(a)-num_unique_els,3):aligns(block_align2(a)-num_unique_els,4));
        
        [block1, block2, ~, ~, ~, ~] = NW_alignment(block1,block2,NW_params);
        align1 = strcat(align1,block1);
        align2 = strcat(align2,block2);
        disp_blocks_1 = strcat(disp_blocks_1,repmat('|',1,length(block1)));
        disp_blocks_2 = strcat(disp_blocks_1,repmat('|',1,length(block2)));
    end
end

disp(disp_blocks_1)
disp(align1)
disp(align2)
disp(disp_blocks_2)

% Add underscores to make the alignments equal in length:
if length(align1) < length(align2)
    align1 = strcat(align1,repmat('_',1,length(align2)-length(align1)));
end
if length(align2) < length(align1)
    align2 = strcat(align2,repmat('_',1,length(align1)-length(align2)));
end


% find the overlaps:
if show_flags
    disp('Computing final outputs')
end
start_overlap = max(min(find(align1 ~= '_')),min(find(align2 ~= '_')));
end_overlap = min(max(find(align1 ~= '_')),max(find(align2 ~= '_')));

if end_overlap < start_overlap
    overlap = [NaN,NaN];
    frac_overlap = 0;
    frac_align = 0;
else
    overlap = [start_overlap, end_overlap];
    frac_overlap = (end_overlap - start_overlap + 1)/length(align1);
    frac_align = sum(align1(start_overlap:end_overlap) == align2(start_overlap:end_overlap))/(end_overlap - start_overlap + 1);
end

if ~isnan(overlap(1)) && ~isnan(overlap(2))
    if overlap(1) > 1
        bh = size(strfind(align1(1:overlap(1)-1),'_'),2);
    else
        bh = 0;
    end
    overlap1 = [overlap(1)-bh,overlap(2)-bh-size(strfind(align1(overlap(1):overlap(2)),'_'),2)];

    if overlap(1) > 1
        bh = size(strfind(align2(1:overlap(1)-1),'_'),2);
    else
        bh = 0;
    end
    overlap2 = [overlap(1)-bh,overlap(2)-bh-size(strfind(align2(overlap(1):overlap(2)),'_'),2)];
else
    overlap1 = [NaN,NaN];
    overlap2 = [NaN,NaN];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [align1, align2] = NW_block_align(seq1,seq2,match_vals,num_unique_els,NW_params)


match_val = NW_params(1);
shift_penalty = NW_params(2);
indel_penalty = NW_params(3);
mismatch_penalty = NW_params(4);


len1 = size(seq1,2);
len2 = size(seq2,2);

matrix = zeros(len1+1,len2+1);
arrows = zeros(len1+1,len2+1);
arrows(1,:) = 2;
arrows(:,1) = 1;
arrows(1,1) = 0;
for r = 2:len1+1
    matrix(r,1) = (r-1)*shift_penalty;
end
for c = 2:len2+1
    matrix(1,c) = (c-1)*shift_penalty;
end
for c = 2:len2+1
    for r = 2:len1+1
        if seq1(r-1)==seq2(c-1)

            if seq1(r-1) <= num_unique_els
                match_v = match_val;
            else
                match_v = match_val*match_vals(2,seq1(r-1)-num_unique_els);
            end
            
            diag = matrix(r-1,c-1)+match_v;
        else
            if seq1(r-1) <= num_unique_els
                mismatch_pen = mismatch_penalty;
            else
                mismatch_pen = mismatch_penalty*match_vals(2,seq1(r-1)-num_unique_els);
            end
            
            diag = matrix(r-1,c-1)+mismatch_pen;
        end
        seq1_indel = matrix(r-1,c)+indel_penalty;
        seq2_indel = matrix(r,c-1)+indel_penalty;
        
        max_val = max(max(seq1_indel, seq2_indel), diag);
        matrix(r,c) = max_val;
        
        if seq1_indel == max_val
            arrows(r,c) = 1;
        elseif seq2_indel == max_val
            arrows(r,c) = 2;
        else 
            arrows(r,c) = 0;
        end
    end
end

max_score = -inf;
for col = 2:len2+1
    if matrix(len1+1,col) >= max_score
        max_score = matrix(len1+1,col);
        c = col;
        r = len1+1;
    end
end
for row = 2:len1+1
    if matrix(row,len2+1) >= max_score
        max_score = matrix(row,len2+1);
        r = row;
        c = len2+1;
    end
end
start_row = r;
start_col = c;

align1 = [];
align2 = [];
done = false;
while ~done
    align = arrows(r,c);
    if align == 1
        align1 = [seq1(r-1),align1];
        align2 = [0,align2];
        r = r - 1;
    elseif align == 2 
        align1 = [0,align1];
        align2 = [seq2(c-1),align2];
        c = c - 1;
    else
        align1 = [seq1(r-1),align1];
        align2 = [seq2(c-1),align2];
        r = r - 1;
        c = c - 1;
    end
    
    if r == 1 && c == 1
        done = true;
    end
end


if start_row <= len1
    align1 = [align1, seq1(start_row:len1)];
    align2 = [align2, repmat(0,1,length(align1)-length(align2))];
end
if start_col <= len2
    align2 = [align2, seq2(start_col:len2)];
    align1 = [align1, repmat(0,1,length(align2)-length(align1))];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [reps,align1,align2] = find_nonoverlapping_alignments(seq1,seq2,min_len_inrep,match_len,indel_len)


% INPUTS: 
% seq1: a [1xm] character vector of nucleotides
% seq2: a [1xn] character vector of nucleotides

% OUTPUTS:
% inreps: a [ix2] cell array, where i is the number of aligned sequences; 
% each of the i rows has format:
%   inreps{i,1}(1:2): the range on seq1 of an aligned sequence
%   inreps{i,2}(1:2): the range on seq2 of corresponding aligned sequence

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
            found = true;

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


for r = 1:size(reps,1)
    % matches:
    fprintf('\n[%d, %d] and [%d, %d]\n',reps{r,1}(1),reps{r,1}(2),reps{r,2}(1),reps{r,2}(2))
    alignseq1 = seq1(reps{r,1}(1):reps{r,1}(2));
    alignseq2 = seq2(reps{r,2}(1):reps{r,2}(2));
    if ~strcmp(alignseq1,alignseq2)
        [alignseq1,alignseq2,~,~,~,~] = NW_alignment(alignseq1,alignseq2,[1,0,0,0]);
    end
    non_matches = 0;
    matchstr = repmat(' ',1,length(align1));
    for m = 1:min(length(alignseq1),length(alignseq2))
        if strcmp(alignseq1(m),'_') || strcmp(alignseq2(m),'_') || ~strcmp(alignseq1(m),alignseq2(m))
            non_matches = non_matches + 1;
            matchstr(m) = 'X';
        end
    end
    disp(alignseq1)
    disp(alignseq2)
    disp(matchstr)
    fprintf('%d\n',non_matches)
end


figure
hold on
axis([0 max(len1,len2) 0 3])

for r = 1:size(reps,1)
    L1 = reps{r,1}(1);
    R1 = reps{r,1}(2);
    midpt1 = L1+(R1-L1)/2;
    L2 = reps{r,2}(1);
    R2 = reps{r,2}(2);
    midpt2 = L2+(R2-L2)/2;
    fill([R1,L1,L1,R1],[1,1,0,0],'b','FaceAlpha',0.2)
    fill([R2,L2,L2,R2],[3,3,2,2],'b','FaceAlpha',0.2)
    line([midpt1,midpt2],[1,2],'Color','b')
end