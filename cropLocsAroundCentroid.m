function cropped_locs = cropLocsAroundCentroid(locs,x_c,y_c,windowsize)

cropped_locs = locs;
cropped_locs = cropped_locs(cropped_locs.x > y_c - (windowsize/2),:);
cropped_locs = cropped_locs(cropped_locs.x < y_c + (windowsize/2),:);
cropped_locs = cropped_locs(cropped_locs.y > x_c - (windowsize/2),:);
cropped_locs = cropped_locs(cropped_locs.y < x_c + (windowsize/2),:);        

end