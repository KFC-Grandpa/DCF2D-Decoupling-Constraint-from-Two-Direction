function [MatingPool] = MatingSelection(Population,type)

    MatingPool = [];

    N = size(Population,2);
    
    PN = N/100;
    if ~mod(PN,1)
        RN = randi(PN)-1;
    else
        RN = 0;
    end



    Density =  DensityCal(Population.objs);
    start = 1+RN * 100;
    
    if (RN+1) * 100 > size(Population,2)
        last = size(Population,2);
    else
        last = (RN+1) * 100;
    end
    index1 = randi([1+RN * 100,last],1,200);
    index2 = randi([1+RN * 100,last],1,200);
    for i = 1:200
        PopObj = Population.objs;
        rand1 = rand;
        rand2 = rand;
        if type == 1 %CDP
            if  PopObj(index1(i),:)<PopObj(index2(i),:)
                MatingPool = [MatingPool, index1(i)];
            elseif  PopObj(index1(i),:)>PopObj(index2(i),:)
                MatingPool = [MatingPool, index2(i)];
            elseif Density(index1(i)) < Density(index2(i))
                MatingPool = [MatingPool, index1(i)];
            else
                MatingPool = [MatingPool, index2(i)];
            end
        elseif  type == 2 %密度优先
            if  Density(index1(i)) < Density(index2(i))
                MatingPool = [MatingPool, index1(i)];
            else
                MatingPool = [MatingPool, index2(i)];
            end
        elseif type == 3 %reverse
            if  PopObj(index1(i),:)<PopObj(index2(i),:)
                MatingPool = [MatingPool, index2(i)];
            elseif  PopObj(index1(i),:)>PopObj(index2(i),:)
                MatingPool = [MatingPool, index1(i)];
            elseif Density(index1(i)) < Density(index2(i))
                MatingPool = [MatingPool, index1(i)];
            else
                MatingPool = [MatingPool, index2(i)];
            end
        end
    end
end




 function Density = DensityCal(PopObj)
    [N,~] = size(PopObj);  
    Zmin       = min(PopObj,[],1);
    PopObj = (PopObj-repmat(Zmin,N,1))./(repmat(max(PopObj),N,1)-repmat(Zmin,N,1)+1e-20);
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = Inf;
    DistanceSort = sort(Distance,2);
    Density = 1./(1+DistanceSort(:,floor(sqrt(length(Distance)))+1));  
 end
