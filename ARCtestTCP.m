%% SERVER SESSION

[~,hostname] = system('hostname');
hostname = string(strtrim(hostname));
address = resolvehost(hostname,"address");

%%

server = tcpserver(address,5000,"ConnectionChangedFcn",@connectionFcn);

%%

configureCallback(server,"byte",7688,@readDataFcn);

%% CLIENT SESSION

client = tcpclient(server.ServerAddress,server.ServerPort,"Timeout",5);
pause(1);

%%

rawData = read(client,961,"double");
reshapedData = reshape(rawData,31,31);
surf(reshapedData);

%%

write(client,rawData,"double");
