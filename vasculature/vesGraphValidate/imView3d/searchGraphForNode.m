function targetFound = searchGraphForNode( E, pt, num, ptTarget)

if num==0
    targetFound = 0;
    return
end

targetFound = 0;

lst = find( E(:,1)==pt | E(:,2)==pt );

for ii=1:length(lst)
    jj = lst(ii);
        
    En = E;
    En(jj,:) = 0;

    v2 = setdiff(E(jj,:),pt);   
    
    if v2==ptTarget | targetFound==1
        targetFound = 1;
        return
    else
        lstE = find(En(:,1)==v2 | En(:,2)==v2);
        nConn1 = setdiff(unique(En(lstE,:)),v2);
    
        targetFound = searchGraphForNode( En, v2, num-1, ptTarget);
    end
    
end

   