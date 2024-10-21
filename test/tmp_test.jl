using Logging

# 打开一个文件用于写入日志
log_file = open("logfile.log", "w")
# 创建一个简单的 logger，将日志写入文件
logger = SimpleLogger(log_file)

# 将全局 logger 设置为我们创建的 logger
global_logger(logger)

for i in 1:10
    @info "This $i message will be written to the log file"
end

# 完成后关闭文件
close(log_file)