function lk=stahlLogLk_quad(events, L, nu, p)
% Calculate log likelihoods according to the Stahl model,
% but assuming that events are observed in a quartet.
%
% Events is a cell array, with one entry per (parental) individual.
% Each entry is a (vertical) vector of event positions (in Morgans).
%
% max map is the chromosome length (in Morgans).
% nu is the interference parameter
% p is the escape probability
%
% Note this function is likely to be very difficult to make sense of 
% as it has been highly optimized.

    disp(['nu = ' num2str(nu) ' p = ' num2str(p)])
    if (nu <= 0)
        lk = -99999999+nu;
        return;
    end
	if (nu > 100)
        lk = -99999999-nu;
        return;
	end
    if ((p < 0) || (p > 1))
        lk = -99999999+p;
        return;
    end
    
    p = max(p, eps);

    [N_indv, N_chr] = size(events);
    assert(N_chr == length(L), 'events doesn''t match L in size.');
    N_tot = N_indv * N_chr;

    two_p = 2*p;
    lk = 0;
    
    log_1_minus_Gstar = zeros(N_tot, 1);
    log_1_minus_Gstar_escape = zeros(N_tot, 1);
    max_map = zeros(N_tot, 1);
    all_events = events(:);
    
    tmp1 = zeros(N_chr, 1);
    tmp2 = zeros(N_chr, 1);
    parfor c=1:N_chr
        tmp1(c) = log(1-Gstar2(L(c), nu, p));
        tmp2(c) = log(1-Gstar1(L(c), two_p));
    end
    
    for c=1:N_chr
        idx = ((c-1)*N_indv)+1:((c-1)*N_indv)+N_indv;
        log_1_minus_Gstar(idx) = tmp1(c);
        log_1_minus_Gstar_escape(idx) = tmp2(c);
        max_map(idx) = L(c);
    end
    
    N_partitions = cellfun(@(x)3^length(x), all_events);
    [uniq_partitions] = unique(N_partitions);
    
    idx = (N_partitions == 1);
    lk = lk + sum(2*log_1_minus_Gstar(idx) + log_1_minus_Gstar_escape(idx));
    
    % Loop over partitions, grouping by the ones with the same number
    parfor u=2:length(uniq_partitions)
        N_part = uniq_partitions(u);    % The number of possible partitions
        idx = (N_partitions == N_part);
        
        partitions = dec2base3(0:(N_part-1));

        event_subset = events(idx);
        if (N_indv == 1)
            all_pos = cell2mat(event_subset');
        else
            all_pos = cell2mat(event_subset);
        end

        lka = repmat(log_1_minus_Gstar(idx), 1, N_part);
        lkc = repmat(log_1_minus_Gstar_escape(idx), 1, N_part);
        max_map_idx = max_map(idx);

        I0 = (partitions == 0);
        I1 = (partitions == 1);
        I2 = (partitions == 2);

        % Only consider unique partitions.
        [uniq_I0,ia0,ib0] = unique(I0, 'rows');
        uniq_I0 = uniq_I0';
        lka_prime = lka(:, ia0);

        sum_sites0 = sum(uniq_I0,1);
        N_sites0 = size(uniq_I0, 1);

        uniq_I2 = fliplr(I2(ia0,:)');
        ia2 = flipud(ia0);
        ib2 = flipud(ib0);
        lkc_prime = lkc(:, ia2);

        % Loop over partitions with the same number of sites
        for n=1:N_sites0
            J0 = (sum_sites0 == n);

            uniq_I00 = uniq_I0(:,J0)';
            uniq_I000 = [true(size(uniq_I00, 1),1) uniq_I00 true(size(uniq_I00, 1),1)];

            all_pos0 = repmat([zeros(size(all_pos, 1),1) all_pos max_map_idx], size(uniq_I000, 1), 1)';
            uniq_I0000 = reshape(repmat(uniq_I000', size(all_pos, 1), 1), size(all_pos0));
            all_pos00 = reshape(all_pos0(uniq_I0000), n+2, []);
            x0 = diff(all_pos00)';

            y0 = log(gstar2(x0(:,1), nu, p) .* (1-Fstar2(x0(:,end), nu, p)));
            lka_prime(:,J0) = reshape(y0, size(all_pos,1), []);

            uniq_I22 = uniq_I2(:,J0)';
            uniq_I222 = [true(size(uniq_I22, 1),1) uniq_I22 true(size(uniq_I22, 1),1)];

            all_pos2 = repmat([zeros(size(all_pos, 1),1) all_pos max_map_idx], size(uniq_I222, 1), 1)';
            uniq_I2222 = reshape(repmat(uniq_I222', size(all_pos, 1), 1), size(all_pos2));
            all_pos22 = reshape(all_pos2(uniq_I2222), n+2, []);
            x2 = diff(all_pos22)';

            y2 = log(gstar1(x2(:,1), two_p) .* (1-Fstar1(x2(:,end), two_p)));
            lkc_prime(:,J0) = reshape(y2, size(all_pos,1), []);

            if (n > 1)
                s = sum(J0);
                lka_prime(:,J0) = lka_prime(:,J0) + reshape(log(prod(fstar2(x0(:,2:end-1), nu, p), 2)), [], s);
                lkc_prime(:,J0) = lkc_prime(:,J0) + reshape(log(prod(fstar1(x2(:,2:end-1), two_p), 2)), [], s);
            end
        end

        % Expand the unique partitions back to the complete matrix
        lka = lka_prime(:, ib0);
        lkc = lkc_prime(:, ib2);

        % No need to calculate lkb, as we can derive it from lka by rearranging...
        [~,idx1] = ismember(I1, I0, 'rows');
        lkb = lka(:, idx1);

        tmp = lka + lkb + lkc;
        ma = max(tmp, [], 2);
        lk = lk + sum(sum(bsxfun(@plus, ma, log(sum(exp(bsxfun(@minus, tmp, ma)),2)))));
    end
            
    lk = lk + log(1-p);
end

function out = fstar1(x, p)
    k_max = 25;
    k_div = 0.5.^(1:k_max)';
    a = repmat(x, [1 1 k_max]);
    b = repmat(reshape(1:k_max, [1 1 k_max]), size(x));
    c = repmat(reshape(k_div, [1 1 k_max]), size(x));
    tmp = gampdf(a, b, 0.5/p) .* c;
    out = sum(tmp, 3);
end

function out = Fstar1(x, p)
    k_max = 25;
    k_div = 0.5.^(1:k_max)';
    a = repmat(x, [1 1 k_max]);
    b = repmat(reshape(1:k_max, [1 1 k_max]), size(x));
    c = repmat(reshape(k_div, [1 1 k_max]), size(x));
    tmp = gamcdf(a, b, 0.5/p) .* c;
    out = sum(tmp, 3);
end

function out = gstar1(x,p)
    out = p*(1-Fstar1(x,p));
end

function out = Gstar1(x,p)
    out = integral(@(X)gstar1(X,p), 0, x);
end

function out = fstar2(x,nu,p)
    k_max = 25;
    k_div = 0.5.^(1:k_max)';
    a = repmat(x, [1 1 k_max]);
    b = repmat(reshape((1:k_max)*nu, [1 1 k_max]), size(x));
    c = repmat(reshape(k_div, [1 1 k_max]), size(x));
    tmp = gampdf(a, b, 0.5/((1-p)*nu)) .* c;
    out = sum(tmp, 3);
end

function out = Fstar2(x,nu,p)
    k_max = 25;
    k_div = 0.5.^(1:k_max)';
    a = repmat(x, [1 1 k_max]);
    b = repmat(reshape((1:k_max)*nu, [1 1 k_max]), size(x));
    c = repmat(reshape(k_div, [1 1 k_max]), size(x));
    tmp = gamcdf(a, b, 0.5/((1-p)*nu)) .* c;
    out = sum(tmp, 3);
end

function out = gstar2(x,nu,p)
    out = (1-p)*(1-Fstar2(x,nu,p));
end

function out = Gstar2(x,nu,p)
    out = integral(@(X)gstar2(X,nu,p), 0, x);
end

function s=dec2base3(d)
    d = d(:);
    n = max(1,round(log2(max(d)+1)/log2(3)));
    while any(3.^n <= d)
        n = n + 1;
    end
    s(:,n) = rem(d,3);
    while any(d) && n >1
        n = n - 1;
        d = floor(d/3);
        s(:,n) = rem(d,3);
    end
end


