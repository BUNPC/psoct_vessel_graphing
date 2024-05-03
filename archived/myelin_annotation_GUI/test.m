figure;im=imagesc(squeeze(img(:,:,10)));
d=drawline(gca,'Label','d1');
d2=drawline(gca,'Label','d2');

l = addlistener(d,'ROIMoved',@(src,evt)roievents(src,evt));
l2 = addlistener(d2,'ROIMoved',@(src,evt)roievents(src,evt));

function roievents(~,evt)
    evname = evt.EventName;
    s = evt.Source.Label;
    if strcmp(evname,'ROIMoved')
        msgbox(s);
    end
end