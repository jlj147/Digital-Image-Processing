%Create figure window with the title Project #2 - Jayce Jones
        movegui(figure('Name','Application Login','NumberTitle','off','ToolBar','None','Position',[0 0 200 100]),'center');
        
        %Create a pushbutton to select and open a new source image as the "current" image
        uicontrol('Style','text','String','Invalid Credentials: Log in as Operator','Position',[80,35,100,40]);
        
        uicontrol('Style','pushbutton','String','OK','Position',[105,25,50,20]);