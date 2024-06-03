addpath 'E:\previous documents\毕业论文\程序及数据\Syn-X联合实验\New C method 3\Project';
Project=@(x) ProjectW(x,Iini,W,340);

FL0TV=@(x)LpFilter(x,{[-1,1],[-1;1]},0,0.005,2);
FL1TV=@(x)LpFilter(x,{[-1,1],[-1;1]},1,0.3,2);
FLpTV=@(x)LpFilter(x,{[-1,1],[-1;1]},0.5,0.03,2);

FL0LA=@(x)LpFilter(x,{[-1,2,-1],[-1;2;-1]},0,0.005,2);
FL1LA=@(x)LpFilter(x,{[-1,2,-1],[-1;2;-1]},1,0.3,2);
FLpLA=@(x)LpFilter(x,{[-1,2,-1],[-1;2;-1]},0.5,0.01,2);

FL0w=@(x)LpFilter(x,{[-1,1],[-1;1],[-1,2,-1],[-1;2;-1]},0,0.003,2);
FL1w=@(x)LpFilter(x,{[-1,1],[-1;1],[-1,2,-1],[-1;2;-1]},1,0.1,2);%0.5
FLpw=@(x)LpFilter(x,{[-1,1],[-1;1],[-1,2,-1],[-1;2;-1]},0.3,0.01,2);


R={FL0TV,FL1TV,FLpTV,FL0LA,FL1LA,FLpLA,FL0w,FL1w,FLpw};
name={'FL0TV','FL1TV','FLpTV','FL0LA','FL1LA','FLpLA','FL0w','FL1w','FLpw'};

for i=1:length(R)
    image = ADMMSolve(Iini,R{i},Project,50,[]);
    save(['result/',name{i},'.mat'],'image');
    imwrite(image,['result/',name{i},'.bmp']);  

end


image=cell(1,12);
image{1}=Syn;
image{2}=Cx;
image{3}=Iini;
 for i=4:12
     data=load(['result/',name{i-3},'.mat']);
     image{i}=data.image;
 end