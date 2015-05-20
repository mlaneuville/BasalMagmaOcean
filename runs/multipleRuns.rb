#!/usr/bin/ruby
include Process

$exe = "a.out"
$maxRuns = 8
$running = 0
$rCount = 0

thickness = [200, 300, 400, 500, 600]
liquidusDrop = [500, 1000]
densityDrop = [500, 1000]

runs = thickness.product(liquidusDrop, densityDrop)

def start(d, l, r)
	cmd = "./" + $exe + " #{d} #{l} #{r} 0"
	pid = fork do
		exec cmd
	end
	puts "(run=%03d) %s started (pid=%05d)" % [$rCount, cmd, pid]
	return pid
end

pids = []

while 1 do
	if (pids.length < $maxRuns)
		x = runs[$rCount][0]
		y = runs[$rCount][1]
		z = runs[$rCount][2]
		pids << start(x, y, z)
		$rCount += 1
	else
		pid = wait()
		puts ""
		puts "%s -- %05d finished" % [Time.now, pid]
		pids.delete_if{|p| p == pid}
	end
end
