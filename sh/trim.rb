#!/usr/bin/env ruby

flag = false
File.open('../csv/record_scaler_beam_on.data','r').each do |l|
  if l.include?('11 -')
    flag = true
  end
  if flag and not l.include?('11 -')
    puts l
  end
end
