function  [w1step, w2step, prstep, w02step, extstep] =newstepsize(w1max,w1min,w2max,w2min,prmax,prmin,w02max,w02min,extmax,extmin,lw1,lw2,lpr,lw02,lext) 


    w1step = (w1max-w1min)/lw1;% w1min=0.3; w1max = 0.6; %w1min must be an integer multiple of w1step.
    w2step = (w2max - w2min)/lw2; %w2min = 3.8; w2max = 4.2;
    prstep = (prmax-prmin)/lpr; %prmin=0; prmax = 0.35;
    w02step = (w02max-w02min)/lw02; %w02min = .85; w02max = .96;
    extstep = (extmax-extmin)/lext; %extmin = 0; extmax = 0;

    if w1step == 0
        w1step = 1;
    end

    if w2step == 0
        w2step = 1;
    end

    if prstep == 0
        prstep = 1;
    end

    if w02step == 0
        w02step = 1;
    end

    if extstep == 0
        extstep = 1;
    end

end

