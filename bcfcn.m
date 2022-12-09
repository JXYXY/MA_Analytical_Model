    function res = bcfcn(xa,xb,C1)
        res = [xa(1)
            xb(2)-C1
            xa(3)
            xb(3) ];
    end