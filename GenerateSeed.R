GenerateSeed <- function( nSeeds , SeedStart )
{
  checkinput = 1
  if( ( SeedStart < -2147483647 ) | ( SeedStart > 2147483647 ) | ( SeedStart%%1 != 0 ) )
  {
    cat(" ",fill=TRUE)
    cat("SeedStart should be an integer between -2147483647 and 2147483647" , fill=TRUE )
    cat(" ",fill=TRUE)
    checkinput = 0
  }
  
  if( checkinput==1 )
  {
    set.seed( SeedStart )
    repeat
    {
      Out = ceiling( ( ( runif( nSeeds ) - .5 ) * 2 ) * 2147483647 )
      if( length(unique(Out)) == nSeeds )
      {
        break
      }
    }
    return( Out )
  }
}