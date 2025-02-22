A = rand(3,3);
logFile = fopen('logfile.txt', 'a'); % Mở file ở chế độ append 'a'
fprintf(logFile, 'Matrix A:\n');
for i = 1:size(A,1)
    fprintf(logFile, '%8.4f ', A(i, :)); % Duyệt từng hàng
    fprintf(logFile, '\n'); % Xuống dòng
end
fclose(logFile);
