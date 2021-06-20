% display impulse response functions


close all

if PLOT_C_AGG
    Output = [  C_agg(2:PLOT_PERIODS)];
        
    Leg =   [   'Agg / Idio cons' ];
end

if PLOT_C_IDIO

    Output = [  Output
                C_idio(2:PLOT_PERIODS)];
    Leg =   [  Leg 
               'Idio cons      '  ];
end

        
if PLOT_CAPITAL
    Output = [  Output
                Z_s(1,2:PLOT_PERIODS)
                Z_s(3,2:PLOT_PERIODS)];
            
      Leg = [ Leg
            'Agg cap        ' 
            'Idio cap       '];
end

if PLOT_LABOUR & sigma_n ~= 0 
    %Output = [  Output
    %            N_agg(2:PLOT_PERIODS) 
    %            W_agg(2:PLOT_PERIODS) 
    %            N_idio_dev(2:PLOT_PERIODS)
    %            W_idio(2:PLOT_PERIODS)];
    Output = [  Output
                Controls(i_l,2:PLOT_PERIODS) 
                Controls(i_w,2:PLOT_PERIODS) 
                Controls(i_lsd,2:PLOT_PERIODS)
                Controls(i_wsd,2:PLOT_PERIODS)];
            
   Leg = [ Leg
            'Agg lab        '
            'Agg wage       '
            'Idio lab dev   ' 
            'Idio wage      '];
end
    
if PLOT_CERTAIN
    Output = [  Output
                C_certain(2:PLOT_PERIODS)];
    Leg = [ Leg
            'Certain equiv c'];

end

if PLOT_CERTAIN_DEV
    Diff = zeros(1,PERIODS);
    Diff = C_agg- C_certain;
    Output = [  Output
                Diff(2:PLOT_PERIODS)];
    Leg = [ Leg
            'Certain deviat '];
    disp('Maximimum deviation from certainty equivalance');   
    disp(max(C_agg(:) - C_certain(:)))
end
    

if PLOT_FULL_INF
    Output = [  Output
                C_full(2:PLOT_PERIODS)];
    Leg = [ Leg
            'Full inf c     '];
end

if PLOT_FIRST_ORDER
    Output = [  Output
                Z_s(6,2:PLOT_PERIODS)
                Z_s(7,2:PLOT_PERIODS)
                Z_s(8,2:PLOT_PERIODS)
                Z_s(9,2:PLOT_PERIODS)];
    Leg = [ Leg
            'Estimate of k  '
            'Estimate of a  '
            'Estimate of ks '
            'Estimate of z  '];
end



if PLOT_HIERACHY ;
    Corr = zeros(1,h);
    for i = 1:h
        Output = [  Output
                    Z_s(i*r+PLOT_HIERACHY,2:PLOT_PERIODS)];
    
        temp=corrcoef(Z_s((i-1)*r+PLOT_HIERACHY,2:PLOT_PERIODS),Z_s(i*r+PLOT_HIERACHY,2:PLOT_PERIODS))   ;
        Corr(i) = temp(1,2);
    end
    Corr'
    
end 

hndl=plot(Output');
legend(hndl,Leg);